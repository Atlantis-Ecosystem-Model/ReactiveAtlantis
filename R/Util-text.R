##' @title Parameter file reader
##' @param text Biological parametar file for Atlatnis
##' @param pattern Text that you are looking
##' @param FG Name of the functional groups
##' @param Vector Logic argument, if the data is on vectors or not
##' @param pprey Logic argument, if the data is a pprey matrix or not
##' @param lineal Logic argument,  if the data is in a vector a lineal
##' @return A matrix with the values from the .prm file
##' @author Javier Porobic
##' @export
text2num <- function(text, pattern, FG = NULL, Vector = FALSE, pprey = FALSE, lineal = FALSE){
    if(!isTRUE(Vector)){
        text <- text[grep(pattern = pattern, text)]
        if(length(text) == 0) warning(paste0('\n\nThere is no ', pattern, ' parameter in your file.'))
        txt  <- gsub(pattern = '[[:space:]]+' ,  '|',  text)
        col1 <- col2 <- vector()
        for( i in 1 : length(txt)){
            tmp     <- unlist(strsplit(txt[i], split = '|', fixed = TRUE))
            if(grepl('#', tmp[1])) next
            tmp2    <- unlist(strsplit(tmp[1], split = '_'))
            if(FG[1] == 'look') {
                col1[i] <- tmp2[1]
            } else {
                id.co   <- which(tmp2 %in% FG)
                if(sum(id.co) == 0) next
                col1[i] <- tmp2[id.co]
            }
            col2[i] <- as.numeric(tmp[2])
        }
        if(is.null(FG)) col1 <- rep('FG', length(col2))
        out.t <- data.frame(FG = col1, Value = col2)
        if(any(is.na(out.t[, 1]))){
            out.t <- out.t[-which(is.na(out.t[, 1])), ]
        }
        return(out.t)
    } else {
        l.pat <- grep(pattern = pattern, text)
        nam   <- gsub(pattern = '[[:space:]]+' ,  '|',  text[l.pat])
        fg    <- vector()
        pos   <- 1
        for( i in 1 : length(nam)){
            tmp     <- unlist(strsplit(nam[i], split = '|', fixed = TRUE))
            if(grepl('#', tmp[1]) || (!grepl('^pPREY', tmp[1]) && pprey  == TRUE)) next
            fg[pos] <- tmp[1]
            if(isTRUE(lineal)){
                t.text <- gsub('+[[:space:]]+', ' ',  text[l.pat[i]])
            } else {
                t.text <- gsub('+[[:space:]]+', ' ',  text[l.pat[i] + 1])
            }
            oldw <- getOption("warn")
            options(warn = -1)
            pp.tmp <- matrix(as.numeric(unlist(strsplit(t.text, split = ' +', fixed = FALSE))), nrow = 1)
            options(warn = oldw)
            if(pos == 1) {
                pp.mat <- pp.tmp
            } else {
                if(ncol(pp.mat) != ncol(pp.tmp)) stop('\nError: The pPrey vector for ', tmp[1], ' has ', ncol(pp.tmp), ' columns and should have ', ncol(pp.mat))
                pp.mat <- rbind(pp.mat, pp.tmp)
            }
            pos    <- pos + 1
        }
        if(all(is.na(pp.mat[, 1]))) pp.mat <- pp.mat[, - 1]
        row.names(pp.mat)                  <- fg
        return(pp.mat)
    }
}


##' @title NCtime calculation
##' @param ncfile Netcdf file
##' @return Vector of dates
##' @author Javier Porobic
##' @export
time_calc <- function(ncfile){
    ## Time Vector
    orign   <- unlist(strsplit(ncdf4::ncatt_get(ncfile, 't')$units, ' ', fixed = TRUE))
    if(orign[1] == 'seconds') {
        Time <- ncdf4::ncvar_get(ncfile, 't') / 86400
    } else {
        Time <- ncdf4::ncvar_get(ncfile, 't')
    }
    Time <- as.Date(Time, origin = orign[3])
    return(Time)
}
