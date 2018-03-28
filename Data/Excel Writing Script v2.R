library(openxlsx)
library(devtools)
Sys.setenv(PATH = paste("C:/Rtools/bin", Sys.getenv("PATH"), sep=";"))
Sys.setenv(BINPREF = "C:/Rtools/mingw_$(WIN)/bin/")

##Dataset 5CSRTT##
summary.list = list('3xTG' = map.list$TG3x, '5xFAD'=map.list$TG5x, 'APPPS1'=map.list$APP)

##Dataset PD ##
summary.list = list('3xTG 4 Monthr' = map.list$TG3x$Month4,'3xTG 7 Monthr' = map.list$TG3x$Month7
                    ,'3xTG 10 Monthr' = map.list$TG3x$Month10,'5xFAD 4 Monthr' = map.list$TG5x$Month4
                    ,'5xFAD 7 Monthr' = map.list$TG5x$Month7,'5xFAD 10 Monthr' = map.list$TG5x$Month10
                    ,'APPPS1 4 Monthr' = map.list$APP$Month4,'APPPS1 7 Monthr' = map.list$APP$Month7
                    ,'APPPS1 10 Monthr' = map.list$APP$Month10)

##Dataset PAL ##
summary.list = list('3xTG 4 Month' = map.list$TG3x$Month4,'3xTG 10 Month' = map.list$TG3x$Month10
                    ,'5xFAD 4 Month' = map.list$TG5x$Month4,'5xFAD 10 Month' = map.list$TG5x$Month10
                    ,'APPPS1 4 Month' = map.list$APP$Month4 ,'APPPS1 10 Month' = map.list$APP$Month10)


## Write Excel Sheet ##
write.xlsx(summary.list,file = 'Weston PD Summary Data.xlsx')

