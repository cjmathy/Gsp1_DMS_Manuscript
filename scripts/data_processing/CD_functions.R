
# Parsing function, which returns a list of two dataframes per experiment: data and metadata
parse_CD_file = function(filepath, return_metadata = FALSE) {
  
    connection = file(filepath, "r")
    metadata <- data.frame(param=character(), value=character(), stringsAsFactors=F) 
    last_param = ''
    while (T) {
        line <- strsplit(readLines(connection, n = 1), split = '\t')
        last_param <- line[[1]][1]
        if ( length(line) == 0 | last_param == 'XYDATA') {
            break
        }
        metadata <- rbind(metadata, data.frame('param'=line[[1]][1],
                                               'value'=line[[1]][2],
                                               stringsAsFactors=F))
    }   
    X <- as.character(subset(metadata, param == 'XUNITS')$value)
    Y <- as.character(subset(metadata, param == 'YUNITS')$value)
    Y2 <- as.character(subset(metadata, param == 'Y2UNITS')$value)  
    data <- data.frame(X=character(), Y=character(), Y2=character(), stringsAsFactors=F) 

    while (T) {
        line <- strsplit(readLines(connection, n = 1), split = '\t')
        if ( length(line) == 0 ) {
          break
        }
        data <- rbind(data, data.frame('X'=as.numeric(line[[1]][1]),
                                       'Y'=as.numeric(line[[1]][2]),
                                       'Y2'=as.numeric(line[[1]][3])))
    }
    close(connection)
    colnames(data) <- c(X, Y, Y2)
    
    data <- data %>% 
        mutate('sample' = get_sample_name(filepath)$name,
               'mutant' = unlist(strsplit(sample, ' '))[1],
               'date' = get_sample_name(filepath)$date,
               'position' = get_sample_name(filepath)$position,
               'filename' = basename(filepath)) %>% 
        rename('CD_mdeg' = `CD[mdeg]`,
               'HT_V' = `HT[V]`)
    if (return_metadata) {
        return(list(data = data, metadata = metadata))
    }

    else {
        return(data)
    }

    
  
  return(list(data = data, metadata = metadata))
}

get_sample_name <- function(path) {
    
    # for filenames of the type path='dirpath/YYYYMMDD_{Name(s)}_ExptType.txt
    f <- tools::file_path_sans_ext(basename(path))
    s <- unlist(strsplit(f, '_'))

    if (length(s) == 4) {
        name <- paste(s[2:(length(s)-1)], collapse=' ')
    } else {
        name <- paste(s[2:(length(s)-2)], collapse=' ')
    } 
    
    date <- as.numeric(s[1])
    
    if (s[2] == 'WT') {
        position <- s[2]
    } else {
        position <- substr(s[2], 1, nchar(s[2]) - 1)
    }


    return(list('name' = name, 'position' = position,'date' = date))
}


# Scale CD data from ellipticity (millidegrees) to molar ellipticity (deg*cm*dmol-1)
# p = pathlength (cm), c = concentration (uM)
mdeg_to_mol_ellip <- function(vec, p, c) {    
    return(100 * (vec/1000) / (p * (c * 1e-6)))
}


plot_scan <- function(expt) {
    ggplot(expt$data, aes_string(x = '`NANOMETERS`', y = '`CD[mdeg]`')) +
        geom_line(size = 0.2) +
        scale_y_continuous(breaks= pretty_breaks()) +
        geom_hline(yintercept = 0, color = 'red', linetype = 'dashed', size = 0.2) + 
        ggtitle(paste(expt$name, expt$date, sep='\n')) +
        xlab('Wavelength (nm)') +
        ylab('Ellipticity (mdeg)') +
        theme_classic() +
        theme(
          text = element_text(family = "Helvetica", size = 6),
          axis.title = element_text(size = 6),
          axis.text = element_text(size = 6),
          axis.ticks = element_line(size = 0.05),
          axis.ticks.length = unit(0.05, 'cm'),
          axis.line = element_line(size = 0.1),
          plot.title = element_text(hjust = 0.5, size = 6)
        )
}


# Plotting functions, for melts and scans
plot_melt <- function(expt) {
    
    sigmoid_model <- fit_sigmoid_to_melt(expt$data)
    
    
    if (is.null(sigmoid_model)) {
        print('fit failed')
        data_to_plot <- mutate(expt$data, 'fit' = mean(expt$data$`CD[mdeg]`))
        Tm = 95
    }
    else {
        data_to_plot <- mutate(expt$data, 'fit' = predict(sigmoid_model, `Temperature [C]`))         
        Tm <- summary(sigmoid_model)$coefficients[2,1]
    }
    
  x <- data_to_plot$`Temperature [C]`
  y <- data_to_plot$fit

  plot <-
    ggplot(data_to_plot) +
    geom_point(aes_string(x = '`Temperature [C]`', y = '`CD[mdeg]`'), size = 0.1) +
    geom_line(aes_string(x = '`Temperature [C]`', y = 'fit'), size = 0.2, color = ucsf_colors$pink1) +
    scale_y_continuous(breaks = pretty_breaks(), label=scientific_format()) +
    geom_vline(xintercept = Tm, color = ucsf_colors$blue1, linetype = 'dashed', size = 0.2) +
    ggtitle(paste(expt$name, expt$date, paste0('approx Tm = ', round(Tm, digits = 1)), sep='\n')) +
    xlab('Temperature (C)') +
    ylab('Ellipticity (mdeg)') +
    theme_classic() +
    theme(
      text = element_text(family = "Helvetica", size = 6),
      axis.title = element_text(size = 6),
      axis.text = element_text(size = 6),
      axis.ticks = element_line(size = 0.05),
      axis.ticks.length = unit(0.05, 'cm'),
      axis.line = element_line(size = 0.1),
      plot.title = element_text(hjust = 0.5, size = 6)
    )
  return(list(plot = plot))
}

# function to estimate the derivate of a function y=f(x) using a sliding window (w)
estimate_derivative <- function(x, y, w) {
    y_prime <- c(rep(NA, length(y)))
    for (i in seq(1+w, length(y)-w)) {
      y_prime[i] <- (y[i+w]-y[i-w])/(x[i+w]-x[i])
    }
    return(y_prime)
}


# Function for fitting a sigmoid to a CD melt curve
fit_sigmoid_melt <- function(df, x, y, id) {
    
    x <- pull(df, !!x)
    y <- pull(df, !!y)
    y1 <- estimate_derivative(x, y, w=15)
    approx_Tm <- x[which(y1 == max(y1, na.rm=TRUE))]

    # for the starting values, m = 2000 is good for WT
    #   - m of 1000-2000 is good for WT like slope. Less cooperative folders will have lower m
    #   - mu and mf are negative for decreasing baselines as temperature increases
    #     at 100% unfolded or 100% folded
    tryCatch({
        nls(y ~ ((yf + mf*x) + (yu + mu*x)*exp(m*(1/Tm - 1/x))) / (1 + exp(m*(1/Tm - 1/x))),
            start = list(m = 2000, Tm = approx_Tm, mf = 0, mu = 0, yf = 0, yu = 1),
            lower = c(m = 0, Tm = 50, mf = -2, mu = -2, yf = -3, yu = -Inf),
            upper = c(m = 3000, Tm = 100, mf = 2, mu = 2, yf = Inf, yu = max_yu),
            algorithm = 'port'
            )
    }, error = function(e) {
        print(paste0(unique(id), ' failed, approx_Tm was ', approx_Tm))
        return(NULL)
    })
}
