#' Calculates Spatial Pattern Analysis usisng a T-square sample procedure.
#' @param x - Distance from the random point to the nearest individual
#' @param y - Distance from individual to its nearest neighbor
#' @return Returns the T-Square Index of Spatial Pattern (C); z-score of C;
#' the Distance Index of Dispersion (I); and z-score of I
#' @examples
#' a=c(7, 19, 11, 18, 12, 27, 23, 27, 12, 8, 2, 4, 10, 18, 19, 8, 3, 9, 4, 5)
#' b=c(8, 6, 6, 13, 16, 11, 18, 8, 7, 7, 3, 7, 32, 22, 22, 12, 17, 18, 11, 10)
#' nna(a,b)
#' @references
#' [1] Cottam, G., & Curtis, J. T. (1956). The use of distance measures in phytosociological sampling. Ecology, 37(3), 451-460. doi:10.2307/1930167
#' [2] Diggle, P. J., Besag, J., & Gleaves, J. T. (1976). Statistical analysis of spatial point patterns by means of distance methods. Biometrics, 659-667.
#' [3] Johnson, R. B., & Zimmer, W. J. (1985). A more powerful test for dispersion using distance measurements. Ecology, 66(5), 1669-1675. doi:10.2307/1938029
#' [4] Lamacraft, R. R., Friedel, M. H., & Chewings, V. H. (1983). Comparison of distance based density estimates for some arid rangeland vegetation. Austral Ecology, 8(2), 181-187. doi:10.1111/j.1442-9993.1983.tb01605.x
#' [5] Ludwig, J. A., & Reynolds, J. F. (1988). Statistical ecology: a primer in methods and computing (Vol. 1). John Wiley & Sons.
#' @importFrom stats pnorm
#' @export
nna <-
  function(x, y){
    x2=x^2
    y2=y^2
    C=x2/(x2+(y2/2))
    zC=(mean(C)-0.5)/(sqrt(1/(12*(length(C)))))
    x22=x2^2
    I=(length(x2)+1)*(sum(x22)/(sum(x2)^2))
    zI=(I-2)/(sqrt(4*(length(x2)-1)/((length(x2)+2)*(length(x2)+3))))
    z_score <- function(z, tailed) {
      if(tailed == 1) return(pnorm(-abs(z)))
      if(tailed == 2) return(2*pnorm(-abs(z)))
      if(tailed != 1 | tailed != 2) return('can only be one or two tailed')
    }
    results=c(mean(C), z_score(zC, 2), I, z_score(zI, 2))
    names(results)=c("C", "p-value (zC)", "I", "p-value (zI)")
    return(results)
  }
