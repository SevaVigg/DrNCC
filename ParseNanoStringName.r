ParseNanoStringName <- function(SourceNameStr){

source("R/getRegExpSubStr.r")



numStr 	<- getRegExpSubStr( SourceNameStr, regexpr( "[0-9][0-9]*.RCC", SourceNameStr) ) 	#strings contain .RCC suffix like 9_[09.RCC]
num	<- as.numeric( substr( numStr, 1, nchar(numStr)-4) )					#cut last four letters like 09

batchStr	<- getRegExpSubStr( SourceNameStr, regexpr( "201[0-9]*", SourceNameStr ) )

#:wtypes <- c("IP", "MC", "m618_sox10", "mitfa w2", "tails", "standard")

if( regexpr( "/(IP|ip)\\.", SourceNameStr) != -1){CellType <- "IP"}else{
 if( regexpr( "/(MC|mc)\\.", SourceNameStr) != -1){CellType <- "MC"}else{
  if( regexpr( "m618", SourceNameStr) != -1){CellType <- "m618"}else{
   if( regexpr( "w2", SourceNameStr) != -1){CellType <- "mitfa w2"}else{
    if( regexpr( "tails", SourceNameStr) != -1){CellType = "tails"}else{ CellType <- "standard"}    
     }}}}

if( batchStr == "20170531"){CellType = "m618"}


hpfStr	<- getRegExpSubStr( SourceNameStr, regexpr( "[0-9][0-9] *[hH][pP][fF]", SourceNameStr ) )  #strings contain hpf suffix
hpf	<- as.numeric( substr( hpfStr, 1, nchar( hpfStr)-3) )
if( CellType == "IP" | CellType == "MC") {hpf <- 48}						   # correction for IP and MC 	


if( regexpr( "[0-3][0-9][-_][0-1][0-9][-_]1[67]", SourceNameStr) != -1) {
	dateEx	<- getRegExpSubStr( SourceNameStr, regexpr("[0-3][0-9][-_][0-1][0-9][-_]1[67]", SourceNameStr))  # like "03[-_]05[-_]16"
	dateEx1 <- getRegExpSubStr( dateEx, regexpr( "^[0-3][0-9]", dateEx ) )					 # like "03"
	dateEx2 <- getRegExpSubStr( dateEx, regexpr( "^[0-3][0-9][-_][0-1][0-9]", dateEx ) )			 # like "03[-_05"
	dateEx2 <- getRegExpSubStr( dateEx, regexpr( "[0-1][0-9]*$", dateEx2 ) )				 # like "05"
	dateEx3 <- getRegExpSubStr( dateEx, regexpr( "[0-1][0-9]*$", dateEx ) )					 # like "16"
}else{	
	dateEx  <- getRegExpSubStr( SourceNameStr, regexpr("[ ][ ]*[0-3][0-9][-_](0[1-9]|1[0-2])([ _]|\\.)[^R]", SourceNameStr) ) # like " 15[-_]05
	dateEx  <- getRegExpSubStr( dateEx, regexpr( "[0-3][0-9]*[-_][0-1][0-9]*", dateEx ) )			 # like "15[-_]05"
	dateEx1 <- getRegExpSubStr( dateEx, regexpr( "^[0-3][0-9]*", dateEx ) )					 # like "15"
	dateEx2 <- getRegExpSubStr( dateEx, regexpr( "[0-1][0-9]*$", dateEx ) )
	dateEx3 <- "2016"
}					 				   					# like "05"
if (dateEx1 == "" & dateEx2 == ""){dateEx <- NA }else{dateEx<- paste( dateEx1, "-", dateEx2, "-", dateEx3, sep = "")}


if( regexpr( "sample[s]*", SourceNameStr) != -1){
	    DeskStr	<- getRegExpSubStr( SourceNameStr, regexpr( "sample[s]*[ ][ ]*[0-9][0-9]*-[0-9][0-9] *", SourceNameStr) )  	# like "sample[s] 1-2"
	    DeskStr	<- getRegExpSubStr( DeskStr, regexpr( "[0-9][0-9]*-[0-9][0-9]*", DeskStr) )					# like 1-12
	    DeskStr	<- paste("samples", DeskStr)
  }else{
   if ( regexpr( "mit[tf]a w2", SourceNameStr) != -1){
	   DeskStr	<- getRegExpSubStr( SourceNameStr, regexpr( "(strip|STRIP) *[a-zA-Z]", SourceNameStr) )			
    }else{ 	
   if( regexpr( "mixed", SourceNameStr) != -1){DeskStr <- "mixed"}else{	
    if ( regexpr( "(m618|sox10|plate [0-9]).*([Ss]trip[s]?[ ][a-zA-Z]|plate[ ]?[0-9])", SourceNameStr) != -1){
      DeskStr <-  getRegExpSubStr( SourceNameStr, regexpr( "(m618|sox10|plate [0-9]).*([Ss]trip[s]?[ ][a-zA-Z]|plate[ ]?[0-9])",  SourceNameStr) )
      }else{
	if( regexpr("[Ss][Tt][Rr][Ii][Pp][s]?[ ]?[a-zA-Z0-9].*tails", SourceNameStr) != -1){
	DeskStr <- getRegExpSubStr( SourceNameStr, regexpr("[Ss][Tt][Rr][Ii][Pp][s]?[ ]?[a-zA-Z0-9].*tails", SourceNameStr) )
	}else{
       	if( regexpr("[Ss][Tt][Rr][Ii][Pp][s]?[ ]?[a-zA-Z0-9]", SourceNameStr) != -1){
	DeskStr <- getRegExpSubStr( SourceNameStr, regexpr( "[Ss][Tt][Rr][Ii][Pp][s]?[ ]?[a-zA-Z0-9]", SourceNameStr) )
       }else{		
    if( regexpr( "plate[ ]?[0-9]?", SourceNameStr) != -1){
     DeskStr	<-  getRegExpSubStr( SourceNameStr, regexpr( "plate[ ]?[0-9]?", SourceNameStr) )
    }else{		
     if( regexpr( " new ", SourceNameStr) != -1){
       DeskStr <- getRegExpSubStr( SourceNameStr, regexpr( "new ", SourceNameStr) ) 		# like new
	}else{
	  if( regexpr( " indexed ", SourceNameStr) != -1){
          DeskStr <- getRegExpSubStr( SourceNameStr, regexpr( "indexed ", SourceNameStr) ) 		# like indexed
	  }else{
	   if( regexpr(	" [0-9][0-9]?[-_][0-9][0-9](_|\\.)", SourceNameStr) != -1){
	     DeskDataStr	<-  getRegExpSubStr( SourceNameStr, regexpr( " [0-9][0-9]?[-_][0-9][0-9](_|\\.)", SourceNameStr) )		#like "1-12. "
	     DeskStr		<-  getRegExpSubStr( DeskDataStr, regexpr( "[0-9][0-9]?[-_][0-9][0-9]?", DeskDataStr) )
	   }else{	 	 	
    	     if( regexpr( "[0-9][0-9]*[-_][0-9][0-9]*[ ][ ]*[0-9][0-9]*[-_][0-9][0-9]*", SourceNameStr) != -1){
	       cat("Warning: the date and the sample description have the same format\n")	
	       DeskDateStr	<- getRegExpSubStr( SourceNameStr, regexpr( "[0-9][0-9]*[-_][0-9][0-9]*[ ][ ]*[0-9][0-9]*[-_][0-9][0-9]*", SourceNameStr) ) 	# like 04-06 1-12
	       DeskStr	<- getRegExpSubStr( DeskDateStr, regexpr( "[0-9][0-9]*[-_][0-9[0-9]*$", DeskDateStr))				# like 1-12
    	       }else{
		 if( regexpr( " study[_ ]?[0-9][_ ]", SourceNameStr) != -1){
		 DeskStr	<- getRegExpSubStr( SourceNameStr, regexpr( " study[_ ]?[0-9][_ ]", SourceNameStr) ) 				# like study 4
		 }else{ DeskStr	<- NA}
		     }		         	    
                 }
              } 
	    }
	  } 
         }
	} 
     }   
   }
  }
}

#Exceptions

if( batchStr == "20161012"){ hpf <- 24}								   # correction for 20161012 _bath 24hpf 22-09-16 strip C.csv with wrong 												    # data in cell records	
if( batchStr == "20160301" | batchStr == "20160216"){ DeskStr <- "flawed"; dateEx <- "00-00-2016"}



ans		<- list()
	ans$num 	<- num
	ans$batch	<- batchStr
	ans$hpf		<- hpf
	ans$dateEx	<- dateEx
	ans$Desk	<- DeskStr
	ans$CellType	<- CellType
  
ans		<- lapply(ans, function(x) {if( !is.na(x) & x == ""){x <- NA}else{x}})				#NAs to all fields that are not assigned

ans

}

