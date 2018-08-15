#' @export
setGeneric("supportCppAccess", function(x) standardGeneric("supportCppAccess"))

#' @export
setMethod("supportCppAccess", "ANY", function(x) FALSE)
