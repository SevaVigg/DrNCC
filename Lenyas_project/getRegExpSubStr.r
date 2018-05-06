getRegExpSubStr <- function(sourceStr, regexprRes){		
# this function takes the string 'sourceStr' and returns the substring identified by regexpRes made by regexpr()

ans <-  substr( sourceStr, regexprRes, regexprRes + attr(regexprRes, "match.length") -1)
ans
}

