#' Use the base R load function but name the object(s) being loaded in
#' 
#' @param file file name
#' 
#' @export
#' @return data loaded in
new_load <- function(file) 
{
	temp_space <- new.env()
	var <- load(file, temp_space)
	out <- mget(var, temp_space)
	if (length(out) == 1) out <- out[[1]]
	rm(temp_space)
	return(out)
}
