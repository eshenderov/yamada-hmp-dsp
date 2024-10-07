# Project:      Routine cold storage leads to hyperacute graft loss in pig-to-primate kidney xenotransplantation; hypothermic machine perfusion may be preferred preservation modality in xenotransplantation
# Description:  Methods for Normalization
# Author:       Adam Luo

setMethod(
  "normalize", "NanoStringGeoMxSet",
  function(object, norm_method = c("q3", "quantile"),
           fromElt = "exprs", toElt,
           housekeepers = HOUSEKEEPERS, ...) {
    norm_method <- match.arg(norm_method)
    switch(norm_method,
           "q3" = {
             q3Norm(object,
                    toElt = toElt, fromElt = fromElt, ...
             )
           },
           "quantile" = {
             quantileNorm(
               object,
               toElt = toElt, fromElt = fromElt, ...
             )
           },
    )
  }
)

q3Norm <- function(object, toElt, fromElt) {
  qs <- apply(exprs(object), 2, function(x) stats::quantile(x, probs = 0.75, type = 6))
  pData(object)[["q3normFactors"]] <- qs / ngeoMean(qs)
  assayDataElement(object, toElt, validate = TRUE) <- sweep(assayDataElement(object, fromElt), 2L, qs / ngeoMean(qs), FUN = "/")
  return(object)
}

quantileNorm <- function(object, toElt, fromElt) {
  assayDataElement(object, toElt, validate = TRUE) <- normalize.quantiles(exprs(object))
  return(object)
}
