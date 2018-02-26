library(openxlsx)
library(devtools)
Sys.setenv(PATH = paste("C:/Rtools/bin", Sys.getenv("PATH"), sep=";"))
Sys.setenv(BINPREF = "C:/Rtools/mingw_$(WIN)/bin/")

##Dataset 5CSRTT & PD ##
summary.list = list('3xTG 4 Month' = map.list$TG3x$Month4$Summary,'3xTG 7 Month' = map.list$TG3x$Month7$Summary
                    ,'3xTG 10 Month' = map.list$TG3x$Month10$Summary,'5xFAD 4 Month' = map.list$TG5x$Month4$Summary
                    ,'5xFAD 7 Month' = map.list$TG5x$Month7$Summary,'5xFAD 10 Month' = map.list$TG5x$Month10$Summary
                    ,'APPPS1 4 Month' = map.list$APP$Month4$Summary,'APPPS1 7 Month' = map.list$APP$Month7$Summary
                    ,'APPPS1 10 Month' = map.list$APP$Month10$Summary)
pval.list = list('3xTG 4 Month' = map.list$TG3x$Month4$PVal,'3xTG 7 Month' = map.list$TG3x$Month7$PVal
                 ,'3xTG 10 Month' = map.list$TG3x$Month10$PVal,'5xFAD 4 Month' = map.list$TG5x$Month4$PVal
                 ,'5xFAD 7 Month' = map.list$TG5x$Month7$PVal,'5xFAD 10 Month' = map.list$TG5x$Month10$PVal
                 ,'APPPS1 4 Month' = map.list$APP$Month4$PVal,'APPPS1 7 Month' = map.list$APP$Month7$PVal
                 ,'APPPS1 10 Month' = map.list$APP$Month10$PVal)
eta.list = list('3xTG 4 Month' = map.list$TG3x$Month4$Eta,'3xTG 7 Month' = map.list$TG3x$Month7$Eta
                ,'3xTG 10 Month' = map.list$TG3x$Month10$Eta,'5xFAD 4 Month' = map.list$TG5x$Month4$Eta
                ,'5xFAD 7 Month' = map.list$TG5x$Month7$Eta,'5xFAD 10 Month' = map.list$TG5x$Month10$Eta
                ,'APPPS1 4 Month' = map.list$APP$Month4$Eta,'APPPS1 7 Month' = map.list$APP$Month7$Eta
                ,'APPPS1 10 Month' = map.list$APP$Month10$Eta)

##Dataset PAL ##
summary.list = list('3xTG 4 Month' = map.list$TG3x$Month4$Summary,'3xTG 10 Month' = map.list$TG3x$Month10$Summary
                    ,'5xFAD 4 Month' = map.list$TG5x$Month4$Summary,'5xFAD 10 Month' = map.list$TG5x$Month10$Summary
                    ,'APPPS1 4 Month' = map.list$APP$Month4$Summary ,'APPPS1 10 Month' = map.list$APP$Month10$Summary)
pval.list = list('3xTG 4 Month' = map.list$TG3x$Month4$Pval,'3xTG 10 Month' = map.list$TG3x$Month10$Pval
                 ,'5xFAD 4 Month' = map.list$TG5x$Month4$Pval,'5xFAD 10 Month' = map.list$TG5x$Month10$Pval
                 ,'APPPS1 4 Month' = map.list$APP$Month4$Pval ,'APPPS1 10 Month' = map.list$APP$Month10$Pval)
eta.list = list('3xTG 4 Month' = map.list$TG3x$Month4$Eta,'3xTG 10 Month' = map.list$TG3x$Month10$Eta
                ,'5xFAD 4 Month' = map.list$TG5x$Month4$Eta,'5xFAD 10 Month' = map.list$TG5x$Month10$Eta
                ,'APPPS1 4 Month' = map.list$APP$Month4$Eta ,'APPPS1 10 Month' = map.list$APP$Month10$Eta)

summary.list = list('PD' = map.list$Summary)
## Write Excel Sheet ##
write.xlsx(summary.list,file = 'Weston PD Summary Data.xlsx')
write.xlsx(pval.list, file= 'Weston PAL Pval Data.xlsx')
write.xlsx(eta.list, file= 'Weston PAL Eta Data.xlsx')
