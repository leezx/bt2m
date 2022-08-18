
# get R2 from regression
# GET EQUATION AND R-SQUARED AS STRING
# SOURCE: https://groups.google.com/forum/#!topic/ggplot2/1TgH-kG5XMA

lm_eqn <- function(df){
  m <- lm(y ~ x, df);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2,
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));
}

lm_eqn_coef_p <- function (df)
{
  #### calculate the slope and p-value for simple linear regression ####
  # use example: project/scPipeline/fancy_analysis/mouse_human_comparison_ctrl.ipynb
  lmp <- function (modelobject) {
    # get p value of a lm model
    if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
    f <- summary(modelobject)$fstatistic
    p <- pf(f[1],f[2],f[3],lower.tail=F)
    attributes(p) <- NULL
    return(p)
  }
  #
  df$y <- df$value
  df$x <- df$percentage
  m = lm(y ~ x, df)
  p.raw <- lmp(m)
  if (is.na(lmp(m)))
    p.raw <- 1
  p <- signif(p.raw, 3)
  #     if (p.raw < 0.05) {
  #         p <- "p-value < 0.05 **"
  #     }
  #     else {
  #         p1 <- format(as.double(p.raw), digits = 2)
  #         p <- paste("p-value =", p1)
  #     }
  eq <- paste("Slope = ", format(as.double(coef(m)[2]), digits = 2),
              "\n", "P-value = ", p, sep = "")
  # as.character(eq)
  return(c(as.double(coef(m)[2]), p.raw, as.character(eq)))
}

#' scale the matrix
#' @param raw.data raw.data
#' @param max.value max.value after scale
#'
#' @return scale.data
#' @export
#' @examples
#' scale.data <- pre.scale.data(logcounts(sce_HSCR))
#'
pre.scale.data <- function(raw.data=logcounts(sce_HSCR), max.value=2) {
  # max.value <- 2
  scale.data <- t(x = scale(x = t(x = as.matrix(x = raw.data)), center = T, scale = T))
  scale.data[is.na(scale.data)] <- 0
  scale.data[scale.data > max.value] <- max.value
  scale.data[scale.data < -max.value] <- -max.value
  return(scale.data)
}

## 5. float matrix to interger matrix
df_float_to_int <- function(counts){
  counts<-apply(counts,2,function(x) {storage.mode(x) <- 'integer'; x})
  return(counts)
}

## 13. read csv in transform
read.tcsv <- function(file, header=TRUE, sep=",", ...) {

  n = max(count.fields(file, sep=sep), na.rm=TRUE)
  x = readLines(file)

  .splitvar = function(x, sep, n) {
    var = unlist(strsplit(x, split=sep))
    length(var) = n
    return(var)
  }

  x = do.call(cbind, lapply(x, .splitvar, sep=sep, n=n))
  x = apply(x, 1, paste, collapse=sep)
  out = read.csv(text=x, sep=sep, header=header, ...)
  return(out)
}


