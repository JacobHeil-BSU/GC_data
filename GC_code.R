#Set up environment ####
##Install and load packages
pkgs <- c("vegan", "this.path", "brms", "ANCOMBC", "phyloseq")

# Load and install required packages
for (i in pkgs) { #Installs packages if not yet installed
  if (!require(i, character.only = TRUE)) install.packages(i)
  remove(i)
  }

### Working Directory 
setwd(this.path::here())

#Load in Data ####
data_16s <- read.csv("../Data/16s_dada2/Data_16s.csv", header = TRUE, row.names = 1, check.names = FALSE) 
data_ITS <- read.csv("../Data/ITS_dada2/Data_ITS.csv", header = TRUE, row.names = 1, check.names = FALSE)
tax_ITS <- read.csv("../Data/ITS_dada2/ITS_taxonomy_new.csv", row.names = 1)
data_culture <- read.csv("../Data/culturing/Data_comm_culture_GC.csv", row.names = 2, header = TRUE, check.names = FALSE)
data_culture <- data_culture[,2:13]
metadata <- read.csv("../metadata/GC_metadata.csv", row.names = 1, header = TRUE, check.names = FALSE)
metadata$Date <-  factor(metadata$Date, levels = unique(metadata$Date)) #change date to a factor
metadata$inoculant <- factor(metadata$inoculant, c("none","sterile water","b.amy","whole community"))
metadata$Watering <- factor(metadata$Watering, c("watered","drought"))
metadata$TD <- paste(metadata$Treatment,metadata$Date_cont, sep = "")
b_tax <- read.csv("../Data/16s_dada2/tax_DC_first_6-28.csv", header = TRUE, row.names = 1, check.names = FALSE)
f_tax <- read.csv("../Data/ITS_dada2/ITS_taxonomy_new.csv", header = TRUE, row.names = 1, check.names = FALSE)

#Assign new names for ASVs
#16S names
data_16s <- data_16s[,reorder(colnames(data_16s), order(colnames(data_16s)))]
#put ASVs in numeric order
names_16s <- paste("16s",sprintf('%0.4d', 1:ncol(data_16s)), sep = "") #generate new ASV names
namingconv_16s <- data.frame("Original" = colnames(data_16s), "New" = names_16s) 
colnames(data_16s) <- names_16s
b_tax <- b_tax[reorder(row.names(b_tax), order(row.names(b_tax))),]
row.names(b_tax) <- names_16s

#row.names(data_tax) == namingconv_ITS$Original
#row.names(data_tax) <- namingconv_ITS$New

#ITS names
data_ITS <- data_ITS[,reorder(colnames(data_ITS), order(colnames(data_ITS)))]
tax_ITS <- tax_ITS[reorder(row.names(tax_ITS), order(row.names(tax_ITS))),]
colnames(data_ITS) == row.names(tax_ITS)
#put ASVs in numeric order
ITS_names <- paste("ITS",sprintf('%0.4d', 1:ncol(data_ITS)), sep = "") #generate new ASV names
namingconv_ITS <- data.frame("Original" = colnames(data_ITS), "New" = ITS_names) 
colnames(data_ITS) <- ITS_names
row.names(tax_ITS) <- ITS_names
f_tax <- f_tax[reorder(row.names(f_tax), order(row.names(f_tax))),]
row.names(f_tax) <- ITS_names

#Calculate other variables ####
#variance from mean photo at each time point
var_photo <- c()
for (i in unique(metadata$Date_cont)){
  for (j in metadata[metadata$Date_cont == i,]$Plant_ID){
    a <- metadata[metadata$Date_cont == i,]
    a <- a[a$Plant_ID == j,]
    var_photo <- append(var_photo, as.numeric(a$Photo - mean(metadata[metadata$Date_cont == i,]$Photo, na.rm = TRUE)))
    remove(a)
    }
}
metadata$var_photo <- var_photo

#Photo change from t0
t_t0 <- c()
for(i in 1:6){
  t_t0 <- append(t_t0, metadata[metadata$Date_cont == i,]$Photo - metadata[metadata$Date_cont == 1,]$Photo)
}
metadata$t_t0 <- t_t0

#mortality, binary
dead <- matrix(nrow = length(unique(metadata$Plant_ID)), ncol = 6)
row.names(dead) <- unique(metadata$Plant_ID)
for (i in unique(metadata$Plant_ID)){
  print(i)
  a <- metadata[metadata$Plant_ID == i,]$Photo - 3
  a[a <= 0] <- 0
  a[a > 0] <- 1
  print(a)
  b <- max(which(a > 0))
  b[b == -Inf] <- 0
  mort <- c(rep(1, b), rep(0, 5 - b), NA)
  print(mort)
  dead[which(unique(metadata$Plant_ID) == i),] <- mort 
}

metadata <- cbind(metadata, "mortality" = as.numeric(c(dead[,1],dead[,2],dead[,3],dead[,4],dead[,5], dead[,6])))

#mortality, percent at date
mort_prop <- tapply(X = metadata$mortality[1:400], INDEX = metadata$TD[1:400], sum) /
  tapply(X = metadata$mortality[1:400], INDEX = metadata$TD[1:400], length)
mort_prop_df <- data.frame("t" = 0:4,"DS" = mort_prop[11:15], "DB" = mort_prop[1:5], "DW" = mort_prop[16:20], "WS" = mort_prop[31:35], "WB" = mort_prop[21:25], "WW" = mort_prop[36:40], row.names = 1:5)
mort_prop_m <- cbind(mort_prop_df[,3], NA, mort_prop_df[,2], mort_prop_df[,4], mort_prop_df[,6], NA,
                       mort_prop_df[,5], mort_prop_df[,7])
mort_prop_m <- rep(mort_prop_m, each = 10)
mort_prop_m <- c(mort_prop_m, rep(NA, 80))
metadata$mort_prop <- mort_prop_m

plot(DS + .025 ~ t, mort_prop_df, type = "l", col = adjustcolor("cyan", alpha = .25), lwd = 3, 
     ylab = "proportion alive", xlab = "time point", ylim = c(0,1.01), bty = "n")
text(x = 3.9, y = .4, "DS")
lines(WS + .015 ~ t, mort_prop_df, type = "l", col = adjustcolor("darkblue", alpha = .25), lwd = 3)
text(x = 4, y = .99, "WS")
lines(DB + .005 ~ t, mort_prop_df, type = "l", col = adjustcolor("magenta", alpha = .25), lwd = 3)
text(x = 4, y = .64, "DB")
lines(WB - .005 ~ t, mort_prop_df, type = "l", col = adjustcolor("red", alpha = .25), lwd = 3)
text(x = 4, y = .87, "WB")
lines(DW - .015 ~ t, mort_prop_df, type = "l", col = adjustcolor("orange", alpha = .25), lwd = 3)
text(x = 4, y = .53, "DW")
lines(WW - .025 ~ t, mort_prop_df, type = "l", col = adjustcolor("gold", alpha = .25), lwd = 3)
text(x = 3, y = .54, "WW")

#Richness and abundance
metadata$rich_cult <- ifelse(row.names(metadata) %in% row.names(data_culture), rowSums(data_culture > 0), NA)
metadata$abund_cult <- ifelse(row.names(metadata) %in% row.names(mort_prop[1:5]), row.names(data_culture), rowSums(data_culture), NA)
metadata$ENS_cult <- ifelse(row.names(metadata) %in% row.names(data_culture), round(exp(diversity(data_culture))), NA)
metadata$rich_16s <- ifelse(row.names(metadata) %in% row.names(data_16s), rowSums(data_16s > 0), NA)
metadata$abund_16s <- ifelse(row.names(metadata) %in% row.names(data_16s), rowSums(data_16s), NA)
metadata$ENS_16s <- ifelse(row.names(metadata) %in% row.names(data_16s), round(exp(diversity(data_16s))), NA)
metadata$rich_ITS <- ifelse(row.names(metadata) %in% row.names(data_ITS), rowSums(data_ITS > 0), NA)
metadata$abund_ITS <- ifelse(row.names(metadata) %in% row.names(data_ITS), rowSums(data_ITS), NA)
metadata$ENS_ITS <- ifelse(row.names(metadata) %in% row.names(data_ITS), round(exp(diversity(data_ITS))), NA)
metadata$water_weight_g <- metadata$Whole_Plant_wet_weightg - metadata$Dry_Weight

#Correlations and PCAs ####
#Time-series response
df_PCA_resp_all <- metadata[,c("Pot_Weight","CA","Photo")]
df_PCA_resp_all <- na.omit(df_PCA_resp_all)
df_PCA_resp_all <- scale(df_PCA_resp_all)
corr_resp_all <- cor(df_PCA_resp_all)
corr_resp_all
heatmap(corr_resp_all, Rowv = NA, Colv = NA, mar = c(7,7), main = "Correlation between time series response variables", cexRow = .75, cexCol = .75, scale = "none")
#This for loop produces pvalues for corr_resp_all
cors_df <- c()
for (i in colnames(df_PCA_resp_all)){
  plist <- c()
  for (j in colnames(df_PCA_resp_all)){
    pv <- signif(cor.test(df_PCA_resp_all[,i], df_PCA_resp_all[,j])$p.value,3)
    plist <- append(plist, pv)
  }
  cors_df <- cbind(cors_df, plist)
  remove(pv, plist)
}
row.names(cors_df) <- colnames(df_PCA_resp_all)
colnames(cors_df) <- colnames(df_PCA_resp_all)
cors_df
heatmap(cors_df, Rowv = NA, Colv = NA, mar = c(7,7), main = "P-values for time series response correlations", cexRow = .75, cexCol = .75, scale = "none")
#PCA of time-series responses
PCA_resp_all <- princomp(corr_resp_all)
summary(PCA_resp_all)
PCA_resp_all
plot(PCA_resp_all$scores, type = "p", xlim = c(-1,1), ylim = c(-1,1), pch = 16)
text(PCA_resp_all$scores, labels = row.names(PCA_resp_all$scores), pos = 3)

#End point corr and PCA
df_PCA_resp_end <- metadata[,c("CA","Root_Length_cm","Shoot_Length_cm","Whole_Plant_wet_weight","root_weight_wetg","stem_wetg","leaves_wet_g","root_dry_g","stem_dryg","leaves_dry_g","Dry_Weight","Delta15N","Delta13C","Weight_percN","Weight_percC")]
df_PCA_resp_end <- df_PCA_resp_end[401:480,]
df_PCA_resp_end <- na.omit(df_PCA_resp_end)
df_PCA_resp_end <- scale(df_PCA_resp_end)
corr_resp_end <- cor(df_PCA_resp_end)
View(corr_resp_end)
heatmap(corr_resp_end, Rowv = NA, Colv = NA, mar = c(7,7), main = "Correlation between time series response variables", cexRow = .75, cexCol = .75, scale = "none")
#This for loop produces pvalues for corr_resp_all
cors_df_end <- c()
for (i in colnames(df_PCA_resp_end)){
  plist <- c()
  for (j in colnames(df_PCA_resp_end)){
    pv <- signif(cor.test(df_PCA_resp_end[,i], df_PCA_resp_end[,j])$p.value,3)
    plist <- append(plist, pv)
  }
  cors_df_end <- cbind(cors_df_end, plist)
  remove(pv, plist)
}
row.names(cors_df_end) <- colnames(df_PCA_resp_end)
colnames(cors_df_end) <- colnames(df_PCA_resp_end)
View(cors_df_end)
heatmap(cors_df_end, Rowv = NA, Colv = NA, mar = c(7,7), main = "P-values for end point response correlations", cexRow = .75, cexCol = .75, scale = "none")
#PCA of time-series responses
PCA_resp_end <- princomp(corr_resp_end)
summary(PCA_resp_end)
PCA_resp_end
plot(PCA_resp_end$scores, type = "n", xlim = c(-1.5,1.5), ylim = c(-1.5,1.5), pch = 16)
text(PCA_resp_end$scores, labels = row.names(PCA_resp_end$scores), pos = 3, cex =.5)
arrows(x0 = rep(0, nrow(PCA_resp_end$scores)), y0 = rep(0, nrow(PCA_resp_end$scores)),
       x1 = PCA_resp_end$scores[,1], y1 = PCA_resp_end$scores[,2], length = .1)

#Mortality ####
#plot(metadata$Photo ~ metadata$Date_cont,type = "n", ylab = "Photosynthesis", xlab = "Time Point", 
#     main = "Photosynthesis by treatment group", xlim = c(1,5))
#abline(a = 5, b = 0, lwd = 3, lty = 2)
#c_pallette <- c("red","blue","green","violet","orange","cyan")
#for (i in unique(metadata$Treatment)[c(1,3:5,7:8)]){
#  plot_data <- metadata[metadata$Treatment == i,]
#  for (j in 1:length(unique(plot_data$Plant_ID))){
#    a = unique(plot_data$Plant_ID)[j]
#    lines(plot_data[plot_data$Plant_ID == a,]$Photo, col = adjustcolor(c_pallette[j], alpha = .25), lwd = 5)
#    points(plot_data[plot_data$Plant_ID == a,]$Photo, col = c_pallette[j], pch = 19) 
#  }  
#}
#photo change plots, raw
#Drought-bacteria photo plot ####
plot(metadata[metadata$Treatment == "WN",]$Photo ~ metadata[metadata$Treatment == "WN",]$Date_cont,
     type = "n", ylab = "Photosynthesis", xlab = "Time Point", main = "Photosynthesis for the drought-bacteria group", 
     xlim = c(1,5), ylim = c(0,20))
abline(a = 5, b = 0, lwd = 3, lty = 2)
lines(metadata[metadata$Plant_ID == "WNA",]$Photo, col = adjustcolor("red", alpha = .25), lwd = 5)
points(metadata[metadata$Plant_ID == "WNA",]$Photo, col = adjustcolor("red", alpha = .5), pch = 19)
lines(metadata[metadata$Plant_ID == "WNB",]$Photo, col = adjustcolor("blue", alpha = .25), lwd = 5)
points(metadata[metadata$Plant_ID == "WNB",]$Photo, col = adjustcolor("blue", alpha = .5), pch = 19)
lines(metadata[metadata$Plant_ID == "WNC",]$Photo, col = adjustcolor("green", alpha = .25), lwd = 5)
points(metadata[metadata$Plant_ID == "WNC",]$Photo, col = adjustcolor("green", alpha = .5), pch = 19)
lines(metadata[metadata$Plant_ID == "WND",]$Photo, col = adjustcolor("orange", alpha = .25), lwd = 5)
points(metadata[metadata$Plant_ID == "WND",]$Photo, col = adjustcolor("orange", alpha = .5), pch = 19)
lines(metadata[metadata$Plant_ID == "WNE",]$Photo, col = adjustcolor("cyan", alpha = .25), lwd = 5)
points(metadata[metadata$Plant_ID == "WNE",]$Photo, col = adjustcolor("cyan", alpha = .5), pch = 19)
lines(metadata[metadata$Plant_ID == "WNF",]$Photo, col = adjustcolor("violet", alpha = .25), lwd = 5)
points(metadata[metadata$Plant_ID == "WNF",]$Photo, col = adjustcolor("violet", alpha = .5), pch = 19)
lines(metadata[metadata$Plant_ID == "WNG",]$Photo, col = adjustcolor("maroon", alpha = .25), lwd = 5)
points(metadata[metadata$Plant_ID == "WNG",]$Photo, col = adjustcolor("maroon", alpha = .5), pch = 19)
lines(metadata[metadata$Plant_ID == "WNH",]$Photo, col = adjustcolor("royalblue", alpha = .25), lwd = 5)
points(metadata[metadata$Plant_ID == "WNH",]$Photo, col = adjustcolor("royalblue", alpha = .5), pch = 19)
lines(metadata[metadata$Plant_ID == "WNI",]$Photo, col = adjustcolor("seagreen", alpha = .25), lwd = 5)
points(metadata[metadata$Plant_ID == "WNI",]$Photo, col = adjustcolor("seagreen", alpha = .5), pch = 19)
lines(metadata[metadata$Plant_ID == "WNJ",]$Photo, col = adjustcolor("orange4", alpha = .25), lwd = 5)
points(metadata[metadata$Plant_ID == "WNJ",]$Photo, col = adjustcolor("orange4", alpha = .5), pch = 19)

#Watered-bacteria photo plot ####
plot(metadata[metadata$Treatment == "WB",]$Photo ~ metadata[metadata$Treatment == "WB",]$Date_cont,
     type = "n", ylab = "Photosynthesis", xlab = "Time Point", main = "Photosynthesis for the watered-bacteria group", 
     xlim = c(1,5), ylim = c(-5,45))
abline(a = 5, b = 0, lwd = 3, lty = 2)
lines(metadata[metadata$Plant_ID == "WBA",]$Photo, col = adjustcolor("red", alpha = .25), lwd = 5)
points(metadata[metadata$Plant_ID == "WBA",]$Photo, col = adjustcolor("red", alpha = .5), pch = 19)
lines(metadata[metadata$Plant_ID == "WBB",]$Photo, col = adjustcolor("blue", alpha = .25), lwd = 5)
points(metadata[metadata$Plant_ID == "WBB",]$Photo, col = adjustcolor("blue", alpha = .5), pch = 19)
lines(metadata[metadata$Plant_ID == "WBC",]$Photo, col = adjustcolor("green", alpha = .25), lwd = 5)
points(metadata[metadata$Plant_ID == "WBC",]$Photo, col = adjustcolor("green", alpha = .5), pch = 19)
lines(metadata[metadata$Plant_ID == "WBD",]$Photo, col = adjustcolor("orange", alpha = .25), lwd = 5)
points(metadata[metadata$Plant_ID == "WBD",]$Photo, col = adjustcolor("orange", alpha = .5), pch = 19)
lines(metadata[metadata$Plant_ID == "WBE",]$Photo, col = adjustcolor("cyan", alpha = .25), lwd = 5)
points(metadata[metadata$Plant_ID == "WBE",]$Photo, col = adjustcolor("cyan", alpha = .5), pch = 19)
lines(metadata[metadata$Plant_ID == "WBF",]$Photo, col = adjustcolor("violet", alpha = .25), lwd = 5)
points(metadata[metadata$Plant_ID == "WBF",]$Photo, col = adjustcolor("violet", alpha = .5), pch = 19)
lines(metadata[metadata$Plant_ID == "WBG",]$Photo, col = adjustcolor("maroon", alpha = .25), lwd = 5)
points(metadata[metadata$Plant_ID == "WBG",]$Photo, col = adjustcolor("maroon", alpha = .5), pch = 19)
lines(metadata[metadata$Plant_ID == "WBH",]$Photo, col = adjustcolor("royalblue", alpha = .25), lwd = 5)
points(metadata[metadata$Plant_ID == "WBH",]$Photo, col = adjustcolor("royalblue", alpha = .5), pch = 19)
lines(metadata[metadata$Plant_ID == "WBI",]$Photo, col = adjustcolor("seagreen", alpha = .25), lwd = 5)
points(metadata[metadata$Plant_ID == "WBI",]$Photo, col = adjustcolor("seagreen", alpha = .5), pch = 19)
lines(metadata[metadata$Plant_ID == "WBJ",]$Photo, col = adjustcolor("orange4", alpha = .25), lwd = 5)
points(metadata[metadata$Plant_ID == "WBJ",]$Photo, col = adjustcolor("orange4", alpha = .5), pch = 19)

#photo with non-inoculated group
plot(metadata[metadata$Treatment == "DN",]$Photo ~ metadata[metadata$Treatment == "DN",]$Date_cont,
     type = "p", ylab = "Photosynthesis", xlab = "Time Point", main = "Photosynthesis for non-inoculated, drought plants", 
     xlim = c(1,5), ylim = c(0,20))

#Drought-sterile photo plot ####
plot(metadata[metadata$Treatment == "DS",]$Photo ~ metadata[metadata$Treatment == "DS",]$Date_cont,
     type = "n", ylab = "Photosynthesis", xlab = "Time Point", main = "Photosynthesis for the drought-sterile group", 
     xlim = c(1,5), ylim = c(-5,45))
abline(a = 5, b = 0, lwd = 3, lty = 2)
lines(metadata[metadata$Plant_ID == "DSA",]$Photo, col = adjustcolor("red", alpha = .25), lwd = 5)
points(metadata[metadata$Plant_ID == "DSA",]$Photo, col = adjustcolor("red", alpha = .5), pch = 19)
lines(metadata[metadata$Plant_ID == "DSB",]$Photo, col = adjustcolor("blue", alpha = .25), lwd = 5)
points(metadata[metadata$Plant_ID == "DSB",]$Photo, col = adjustcolor("blue", alpha = .5), pch = 19)
lines(metadata[metadata$Plant_ID == "DSC",]$Photo, col = adjustcolor("green", alpha = .25), lwd = 5)
points(metadata[metadata$Plant_ID == "DSC",]$Photo, col = adjustcolor("green", alpha = .5), pch = 19)
lines(metadata[metadata$Plant_ID == "DSD",]$Photo, col = adjustcolor("orange", alpha = .25), lwd = 5)
points(metadata[metadata$Plant_ID == "DSD",]$Photo, col = adjustcolor("orange", alpha = .5), pch = 19)
lines(metadata[metadata$Plant_ID == "DSE",]$Photo, col = adjustcolor("cyan", alpha = .25), lwd = 5)
points(metadata[metadata$Plant_ID == "DSE",]$Photo, col = adjustcolor("cyan", alpha = .5), pch = 19)
lines(metadata[metadata$Plant_ID == "DSF",]$Photo, col = adjustcolor("violet", alpha = .25), lwd = 5)
points(metadata[metadata$Plant_ID == "DSF",]$Photo, col = adjustcolor("violet", alpha = .5), pch = 19)
lines(metadata[metadata$Plant_ID == "DSG",]$Photo, col = adjustcolor("maroon", alpha = .25), lwd = 5)
points(metadata[metadata$Plant_ID == "DSG",]$Photo, col = adjustcolor("maroon", alpha = .5), pch = 19)
lines(metadata[metadata$Plant_ID == "DSH",]$Photo, col = adjustcolor("royalblue", alpha = .25), lwd = 5)
points(metadata[metadata$Plant_ID == "DSH",]$Photo, col = adjustcolor("royalblue", alpha = .5), pch = 19)
lines(metadata[metadata$Plant_ID == "DSI",]$Photo, col = adjustcolor("seagreen", alpha = .25), lwd = 5)
points(metadata[metadata$Plant_ID == "DSI",]$Photo, col = adjustcolor("seagreen", alpha = .5), pch = 19)
lines(metadata[metadata$Plant_ID == "DSJ",]$Photo, col = adjustcolor("orange4", alpha = .25), lwd = 5)
points(metadata[metadata$Plant_ID == "DSJ",]$Photo, col = adjustcolor("orange4", alpha = .5), pch = 19)

#Watered-sterile photo plot ####
plot(metadata[metadata$Treatment == "WS",]$Photo ~ metadata[metadata$Treatment == "WS",]$Date_cont,
     type = "n", ylab = "Photosynthesis", xlab = "Time Point", main = "Photosynthesis for the watered-sterile group", 
     xlim = c(1,5), ylim = c(-5,45))
abline(a = 5, b = 0, lwd = 3, lty = 2)
lines(metadata[metadata$Plant_ID == "WSA",]$Photo, col = adjustcolor("red", alpha = .25), lwd = 5)
points(metadata[metadata$Plant_ID == "WSA",]$Photo, col = adjustcolor("red", alpha = .5), pch = 19)
lines(metadata[metadata$Plant_ID == "WSB",]$Photo, col = adjustcolor("blue", alpha = .25), lwd = 5)
points(metadata[metadata$Plant_ID == "WSB",]$Photo, col = adjustcolor("blue", alpha = .5), pch = 19)
lines(metadata[metadata$Plant_ID == "WSC",]$Photo, col = adjustcolor("green", alpha = .25), lwd = 5)
points(metadata[metadata$Plant_ID == "WSC",]$Photo, col = adjustcolor("green", alpha = .5), pch = 19)
lines(metadata[metadata$Plant_ID == "WSD",]$Photo, col = adjustcolor("orange", alpha = .25), lwd = 5)
points(metadata[metadata$Plant_ID == "WSD",]$Photo, col = adjustcolor("orange", alpha = .5), pch = 19)
lines(metadata[metadata$Plant_ID == "WSE",]$Photo, col = adjustcolor("cyan", alpha = .25), lwd = 5)
points(metadata[metadata$Plant_ID == "WSE",]$Photo, col = adjustcolor("cyan", alpha = .5), pch = 19)
lines(metadata[metadata$Plant_ID == "WSF",]$Photo, col = adjustcolor("violet", alpha = .25), lwd = 5)
points(metadata[metadata$Plant_ID == "WSF",]$Photo, col = adjustcolor("violet", alpha = .5), pch = 19)
lines(metadata[metadata$Plant_ID == "WSG",]$Photo, col = adjustcolor("maroon", alpha = .25), lwd = 5)
points(metadata[metadata$Plant_ID == "WSG",]$Photo, col = adjustcolor("maroon", alpha = .5), pch = 19)
lines(metadata[metadata$Plant_ID == "WSH",]$Photo, col = adjustcolor("royalblue", alpha = .25), lwd = 5)
points(metadata[metadata$Plant_ID == "WSH",]$Photo, col = adjustcolor("royalblue", alpha = .5), pch = 19)
lines(metadata[metadata$Plant_ID == "WSI",]$Photo, col = adjustcolor("seagreen", alpha = .25), lwd = 5)
points(metadata[metadata$Plant_ID == "WSI",]$Photo, col = adjustcolor("seagreen", alpha = .5), pch = 19)
lines(metadata[metadata$Plant_ID == "WSJ",]$Photo, col = adjustcolor("orange4", alpha = .25), lwd = 5)
points(metadata[metadata$Plant_ID == "WSJ",]$Photo, col = adjustcolor("orange4", alpha = .5), pch = 19)

#Drought-wholecomm photo plot ####
plot(metadata[metadata$Treatment == "DW",]$Photo ~ metadata[metadata$Treatment == "DW",]$Date_cont,
     type = "n", ylab = "Photosynthesis", xlab = "Time Point", main = "Photosynthesis for the drought-wholecomm group", 
     xlim = c(1,5), ylim = c(-5,45))
abline(a = 5, b = 0, lwd = 3, lty = 2)
lines(metadata[metadata$Plant_ID == "DWA",]$Photo, col = adjustcolor("red", alpha = .25), lwd = 5)
points(metadata[metadata$Plant_ID == "DWA",]$Photo, col = adjustcolor("red", alpha = .5), pch = 19)
lines(metadata[metadata$Plant_ID == "DWB",]$Photo, col = adjustcolor("blue", alpha = .25), lwd = 5)
points(metadata[metadata$Plant_ID == "DWB",]$Photo, col = adjustcolor("blue", alpha = .5), pch = 19)
lines(metadata[metadata$Plant_ID == "DWC",]$Photo, col = adjustcolor("green", alpha = .25), lwd = 5)
points(metadata[metadata$Plant_ID == "DWC",]$Photo, col = adjustcolor("green", alpha = .5), pch = 19)
lines(metadata[metadata$Plant_ID == "DWD",]$Photo, col = adjustcolor("orange", alpha = .25), lwd = 5)
points(metadata[metadata$Plant_ID == "DWD",]$Photo, col = adjustcolor("orange", alpha = .5), pch = 19)
lines(metadata[metadata$Plant_ID == "DWE",]$Photo, col = adjustcolor("cyan", alpha = .25), lwd = 5)
points(metadata[metadata$Plant_ID == "DWE",]$Photo, col = adjustcolor("cyan", alpha = .5), pch = 19)
lines(metadata[metadata$Plant_ID == "DWF",]$Photo, col = adjustcolor("violet", alpha = .25), lwd = 5)
points(metadata[metadata$Plant_ID == "DWF",]$Photo, col = adjustcolor("violet", alpha = .5), pch = 19)
lines(metadata[metadata$Plant_ID == "DWG",]$Photo, col = adjustcolor("maroon", alpha = .25), lwd = 5)
points(metadata[metadata$Plant_ID == "DWG",]$Photo, col = adjustcolor("maroon", alpha = .5), pch = 19)
lines(metadata[metadata$Plant_ID == "DWH",]$Photo, col = adjustcolor("royalblue", alpha = .25), lwd = 5)
points(metadata[metadata$Plant_ID == "DWH",]$Photo, col = adjustcolor("royalblue", alpha = .5), pch = 19)
lines(metadata[metadata$Plant_ID == "DWI",]$Photo, col = adjustcolor("seagreen", alpha = .25), lwd = 5)
points(metadata[metadata$Plant_ID == "DWI",]$Photo, col = adjustcolor("seagreen", alpha = .5), pch = 19)
lines(metadata[metadata$Plant_ID == "DWJ",]$Photo, col = adjustcolor("orange4", alpha = .25), lwd = 5)
points(metadata[metadata$Plant_ID == "DWJ",]$Photo, col = adjustcolor("orange4", alpha = .5), pch = 19)

#Watered-wholecomm photo plot ####
plot(metadata[metadata$Treatment == "WW",]$Photo ~ metadata[metadata$Treatment == "WW",]$Date_cont,
     type = "n", ylab = "Photosynthesis", xlab = "Time Point", main = "Photosynthesis for the watered-wholecomm group", 
     xlim = c(1,5), ylim = c(-5,45))
abline(a = 5, b = 0, lwd = 3, lty = 2)
lines(metadata[metadata$Plant_ID == "WWA",]$Photo, col = adjustcolor("red", alpha = .25), lwd = 5)
points(metadata[metadata$Plant_ID == "WWA",]$Photo, col = adjustcolor("red", alpha = .5), pch = 19)
lines(metadata[metadata$Plant_ID == "WWB",]$Photo, col = adjustcolor("blue", alpha = .25), lwd = 5)
points(metadata[metadata$Plant_ID == "WWB",]$Photo, col = adjustcolor("blue", alpha = .5), pch = 19)
lines(metadata[metadata$Plant_ID == "WWC",]$Photo, col = adjustcolor("green", alpha = .25), lwd = 5)
points(metadata[metadata$Plant_ID == "WWC",]$Photo, col = adjustcolor("green", alpha = .5), pch = 19)
lines(metadata[metadata$Plant_ID == "WWD",]$Photo, col = adjustcolor("orange", alpha = .25), lwd = 5)
points(metadata[metadata$Plant_ID == "WWD",]$Photo, col = adjustcolor("orange", alpha = .5), pch = 19)
lines(metadata[metadata$Plant_ID == "WWE",]$Photo, col = adjustcolor("cyan", alpha = .25), lwd = 5)
points(metadata[metadata$Plant_ID == "WWE",]$Photo, col = adjustcolor("cyan", alpha = .5), pch = 19)
lines(metadata[metadata$Plant_ID == "WWF",]$Photo, col = adjustcolor("violet", alpha = .25), lwd = 5)
points(metadata[metadata$Plant_ID == "WWF",]$Photo, col = adjustcolor("violet", alpha = .5), pch = 19)
lines(metadata[metadata$Plant_ID == "WWG",]$Photo, col = adjustcolor("maroon", alpha = .25), lwd = 5)
points(metadata[metadata$Plant_ID == "WWG",]$Photo, col = adjustcolor("maroon", alpha = .5), pch = 19)
lines(metadata[metadata$Plant_ID == "WWH",]$Photo, col = adjustcolor("royalblue", alpha = .25), lwd = 5)
points(metadata[metadata$Plant_ID == "WWH",]$Photo, col = adjustcolor("royalblue", alpha = .5), pch = 19)
lines(metadata[metadata$Plant_ID == "WWI",]$Photo, col = adjustcolor("seagreen", alpha = .25), lwd = 5)
points(metadata[metadata$Plant_ID == "WWI",]$Photo, col = adjustcolor("seagreen", alpha = .5), pch = 19)
lines(metadata[metadata$Plant_ID == "WWJ",]$Photo, col = adjustcolor("orange4", alpha = .25), lwd = 5)
points(metadata[metadata$Plant_ID == "WWJ",]$Photo, col = adjustcolor("orange4", alpha = .5), pch = 19)

#photo all data plot
colors <- palette.colors(palette = "Okabe-Ito")[c(1:4,6:9)]
p_meta <- metadata[metadata$Treatment %in% c("DB","DS","DW","WB","WS","WW"),]

plot(Photo ~ Date_cont, data = p_meta, type = "n", xlim = c(1,5))
#rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "grey")
for(i in 1:6){
stripchart(Photo ~ Date_cont, data = p_meta[p_meta$Treatment == unique(p_meta$Treatment)[i],], 
           vertical = TRUE, na.action = na.omit, col = adjustcolor(colors[i], alpha = 1),
           pch = i, add = TRUE, method = "jitter")
plot_model <- lm(Photo ~ Date_cont, data = p_meta[p_meta$Treatment == unique(p_meta$Treatment)[i],])
pm_sum <- summary(plot_model)
error <- c(.5 * pm_sum$coefficients[2] + pm_sum$coefficients[1] - pm_sum$coefficients[3],
           5.5 * pm_sum$coefficients[2] + pm_sum$coefficients[1] - pm_sum$coefficients[3],
           5.5 * pm_sum$coefficients[2] + pm_sum$coefficients[1] + pm_sum$coefficients[3],
           .5 * pm_sum$coefficients[2] + pm_sum$coefficients[1] + pm_sum$coefficients[3])
abline(plot_model, col = adjustcolor(colors[i], alpha = .5), lwd = 3)
polygon(y = error, x = c(.5, 5.5, 5.5, .5), col = adjustcolor(colors[i], alpha = .1),
        border = adjustcolor(colors[i], alpha = .25))
}

plot(CA ~ Date_cont, data = metadata, type = "n", xlim = c(1,5))
#rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "grey")
for(i in 1:8){
  stripchart(CA ~ Date_cont, data = metadata[metadata$Treatment == unique(metadata$Treatment)[i],], 
             vertical = TRUE, na.action = na.omit, col = adjustcolor(colors[i], alpha = 1),
             pch = i, add = TRUE, method = "jitter")
  plot_model <- lm(CA ~ Date_cont, data = metadata[metadata$Treatment == unique(metadata$Treatment)[i],], na.action = na.exclude)
  pm_sum <- summary(plot_model)
  error <- c(.5 * pm_sum$coefficients[2] + pm_sum$coefficients[1] - pm_sum$coefficients[3],
             5.5 * pm_sum$coefficients[2] + pm_sum$coefficients[1] - pm_sum$coefficients[3],
             5.5 * pm_sum$coefficients[2] + pm_sum$coefficients[1] + pm_sum$coefficients[3],
             .5 * pm_sum$coefficients[2] + pm_sum$coefficients[1] + pm_sum$coefficients[3])
  abline(plot_model, col = adjustcolor(colors[i], alpha = .5), lwd = 3)
  polygon(y = error, x = c(.5, 5.5, 5.5, .5), col = adjustcolor(colors[i], alpha = .1),
          border = adjustcolor(colors[i], alpha = .25))
}

ggplot(data = p_meta[1:300,], aes(x = as.factor(Date_cont), y = Photo, color = Treatment, shape = Treatment)) +
  geom_jitter(aes(fill = Treatment), width = .25) +
  geom_smooth(aes(x = Date_cont ,fill = Treatment), method = "lm", alpha = .1) +
  theme_classic()
  
ggplot(data = p_meta[1:300,], aes(x = as.factor(Date_cont), y = CA, color = Treatment, shape = Treatment)) +
  geom_jitter(aes(fill = Treatment), width = .25) +
  geom_smooth(aes(x = Date_cont ,fill = Treatment), method = "lm", alpha = .1) +
  theme_classic()



test(1:5)  

apply(1:5, FUN = test())

#photo change from t0 plots (t0, t1-t0, t2-t0, etc)
# all plots loop ####
#Drought-bacteria photo plot ####
for (i in unique(metadata$Treatment)){
plot(metadata[metadata$Treatment == i,]$t_t0 ~ metadata[metadata$Treatment == i,]$Date_cont,
     type = "n", ylab = "Photo at tn - t1", xlab = "Time Point", 
     main = paste("Change in Photosynthesis for", i, sep = " "), xlim = c(1,5), ylim = c(-30,35))
abline(a = 0, b = 0, lwd = 3, lty = 2)
lines(metadata[metadata$Plant_ID == paste(i, "A", sep = ""),]$t_t0, col = adjustcolor("red", alpha = .25), lwd = 5)
points(metadata[metadata$Plant_ID == paste(i, "A", sep = ""),]$t_t0, col = adjustcolor("red", alpha = .5), pch = 19)
lines(metadata[metadata$Plant_ID == paste(i, "B", sep = ""),]$t_t0, col = adjustcolor("blue", alpha = .25), lwd = 5)
points(metadata[metadata$Plant_ID == paste(i, "B", sep = ""),]$t_t0, col = adjustcolor("blue", alpha = .5), pch = 19)
lines(metadata[metadata$Plant_ID == paste(i, "C", sep = ""),]$t_t0, col = adjustcolor("green", alpha = .25), lwd = 5)
points(metadata[metadata$Plant_ID == paste(i, "C", sep = ""),]$t_t0, col = adjustcolor("green", alpha = .5), pch = 19)
lines(metadata[metadata$Plant_ID == paste(i, "D", sep = ""),]$t_t0, col = adjustcolor("orange", alpha = .25), lwd = 5)
points(metadata[metadata$Plant_ID == paste(i, "D", sep = ""),]$t_t0, col = adjustcolor("orange", alpha = .5), pch = 19)
lines(metadata[metadata$Plant_ID == paste(i, "E", sep = ""),]$t_t0, col = adjustcolor("cyan", alpha = .25), lwd = 5)
points(metadata[metadata$Plant_ID == paste(i, "E", sep = ""),]$t_t0, col = adjustcolor("cyan", alpha = .5), pch = 19)
lines(metadata[metadata$Plant_ID == paste(i, "F", sep = ""),]$t_t0, col = adjustcolor("violet", alpha = .25), lwd = 5)
points(metadata[metadata$Plant_ID == paste(i, "F", sep = ""),]$t_t0, col = adjustcolor("violet", alpha = .5), pch = 19)
lines(metadata[metadata$Plant_ID == paste(i, "G", sep = ""),]$t_t0, col = adjustcolor("maroon", alpha = .25), lwd = 5)
points(metadata[metadata$Plant_ID == paste(i, "G", sep = ""),]$t_t0, col = adjustcolor("maroon", alpha = .5), pch = 19)
lines(metadata[metadata$Plant_ID == paste(i, "H", sep = ""),]$t_t0, col = adjustcolor("royalblue", alpha = .25), lwd = 5)
points(metadata[metadata$Plant_ID == paste(i, "H", sep = ""),]$t_t0, col = adjustcolor("royalblue", alpha = .5), pch = 19)
lines(metadata[metadata$Plant_ID == paste(i, "I", sep = ""),]$t_t0, col = adjustcolor("seagreen", alpha = .25), lwd = 5)
points(metadata[metadata$Plant_ID == paste(i, "I", sep = ""),]$t_t0, col = adjustcolor("seagreen", alpha = .5), pch = 19)
lines(metadata[metadata$Plant_ID == paste(i, "J", sep = ""),]$t_t0, col = adjustcolor("orange4", alpha = .25), lwd = 5)
points(metadata[metadata$Plant_ID == paste(i, "J", sep = ""),]$t_t0, col = adjustcolor("orange4", alpha = .5), pch = 19)
}

#boxplots for end point variables ####
responses <- c("CA","Root_Length_cm","Shoot_Length_cm","Whole_Plant_wet_weightg","root_weight_wetg","stem_wetg","leaves_wet_g","root_dry_g","stem_dryg","leaves_dry_g","Dry_Weight","Delta15N","Delta13C","Weight_percN","Weight_percC")
for (i in responses){
  boxplot(metadata[401:480,i] ~ metadata[401:480,"Treatment"], na.action = na.omit, ylab = i, xlab = "treatment group")
  stripchart(metadata[401:480,i] ~ metadata[401:480,"Treatment"], add = TRUE, vertical = TRUE, pch = 1, 
             method = "jitter")
  }

#Models ####
#time series
metadata_ts <- metadata[metadata$inoculant %in% c("sterile water","b.amy","whole community"),]
metadata_ts$inoculant <- factor(metadata_ts$inoculant, c("sterile water" ,"b.amy", "whole community"))
metadata_ts$Watering <- factor(metadata_ts$Watering, c("watered","drought"))
metadata_ts[metadata_ts$mort_prop %in% 1,]$mort_prop <- .999

#mortality proportion
model_mort <- brm(mort_prop ~  Watering * inoculant + (1|Plant_ID) + (1|Date), 
                   data = metadata_ts, family = "beta", iter = 3000, 
                   control = list(max_treedepth = 10))
summary(model_mort)
mcmc_plot(model_mort)

#mortality and nutrient content
meta_mort <- metadata_ts[301:360,]
meta_mort$mort_prop1 <- metadata_ts$mort_prop[241:300]
meta_mort$Photo <- metadata_ts$Photo[241:300]
C13_mort <- brm(Delta13C ~ mort_prop1 + Watering * inoculant, data = meta_mort, 
                family = "gaussian", iter = 20000, control = list(max_treedepth = 12))
C13_mort_bayes <- C13_mort
C13_mort <- glm(Delta13C ~ mort_prop + Watering * inoculant, data = meta_mort, 
                family = "gaussian")
mcmc_plot(C13_mort)
plot(conditional_effects(C13_mort))
summary(C13_mort)
plot(metadata_ts$Weight_percC[301:360] ~ metadata_ts$mort_prop[241:300])

C13_mort <- brm(Delta15N ~ mort_prop1 + Watering * inoculant, data = eta_mort, 
                family = "gaussian", iter = 3000, control = list(max_treedepth = 10))

#photosynthesis
#set.seed(33)
model_photo <- brm(Photo ~  0 + inoculant * Watering + (1|Plant_ID) + (1|Date), 
                        data = metadata_ts, family = "Gaussian", iter = 3000, 
                        control = list(max_treedepth = 10))
#saveRDS(model_photo, "model_photo.rds")
model_photo <- readRDS("model_photo.rds")
summary(model_photo)
mcmc_plot(model_photo)
plot(model_photo)
plot(conditional_effects(model_photo))
ce_mp <- conditional_effects(model_photo)
plot(ce_mp$Watering$estimate ~ ce_mp$Watering$Watering, ylim = c(10,22))
arrows(x0 = c(1,2), 
       y0 = c(ce_mp$Watering$estimate[1] + ce_mp$Watering$se__[1],
              ce_mp$Watering$estimate[2] + ce_mp$Watering$se__[2]),
       x1 = c(1,2), 
       y1 = c(ce_mp$Watering$estimate[1] - ce_mp$Watering$se__[1],
              ce_mp$Watering$estimate[2] - ce_mp$Watering$se__[2]),
       angle = 90, code = 3)
plot(ce_mp$inoculant$estimate ~ ce_mp$inoculant$inoculant, ylim = c(5,22))
arrows(x0 = c(1,2,3), 
       y0 = c(ce_mp$inoculant$estimate[1] + ce_mp$inoculant$se__[1],
              ce_mp$inoculant$estimate[2] + ce_mp$inoculant$se__[2],
              ce_mp$inoculant$estimate[3] + ce_mp$inoculant$se__[3]),
       x1 = c(1,2,3), 
       y1 = c(ce_mp$inoculant$estimate[1] - ce_mp$inoculant$se__[1],
              ce_mp$inoculant$estimate[2] - ce_mp$inoculant$se__[2],
              ce_mp$inoculant$estimate[3] - ce_mp$inoculant$se__[3]),
       angle = 90, code = 3)

post_photo <- as.data.frame(model_photo) # extract posterior samples from model fit
hist(post_photo$b_inoculantb.amy)
for(i in 2:ncol(post_photo)){
  print(colnames(post_photo[i]))
  print(sum(post_photo[,i] < 0) / length(post_photo[,i]))
  print(sum(post_photo[,i] > 0) / length(post_photo[,i]))
}

#Canopy area
#set.seed(33)
#model_canopy <- brm(CA ~ Watering * inoculant + (1|Plant_ID), 
#                       data = metadata, family = "Gaussian", iter = 4000, 
#                       control = list(max_treedepth = 10))
#saveRDS(model_canopy, "model_canopy.rds")
model_canopy <- readRDS("model_canopy.rds")
summary(model_canopy)
mcmc_plot(model_canopy)
plot(conditional_effects(model_canopy))
post_canopy <- as.data.frame(model_canopy) # extract posterior samples from model fit
hist(post_canopy$b_inoculantb.amy)
for(i in 2:ncol(post_canopy)){
  print(colnames(post_canopy[i]))
  print(sum(post_canopy[,i] < 0) / length(post_canopy[,i]))
  print(sum(post_canopy[,i] > 0) / length(post_canopy[,i]))
}

#end point models
metadata_ep <- metadata[401:480,]
metadata_ph <- metadata[metadata$Date_cont == 5,] #subset live plants at end point
metadata_ph <- metadata_ph[metadata_ph$Photo > 2 | is.na(metadata_ph$Photo),]
metadata_ep <- metadata_ep[metadata_ep$Plant_ID %in% metadata_ph$Plant_ID,]

#set.seed(33)
#model_ca <- brm(CA ~ Watering * inoculant, 
#                       data = metadata_ep, family = "Gaussian", iter = 2000, 
#                       control = list(max_treedepth = 10))
#saveRDS(model_ca, "model_ca.rds")
model_ca <- readRDS("model_ca.rds")
summary(model_ca)
mcmc_plot(model_ca)
plot(conditional_effects(model_ca))

#root length - model not converging, needs work
#set.seed(33)
#model_rl <- brm(Root_Length_cm ~ Watering * inoculant, 
#                       data = metadata_ep, family = "Gaussian", iter = 2000, 
#                       control = list(max_treedepth = 10))
#saveRDS(model_rl, "model_rl.rds")
model_rl <- readRDS("model_rl.rds")
summary(model_rl)
mcmc_plot(model_rl)
plot(conditional_effects(model_rl))
post_rl <- as.data.frame(model_rl) # extract posterior samples from model fit
for(i in 2:ncol(post_rl)){
  print(colnames(post_rl[i]))
  print(sum(post_rl[,i] < 0) / length(post_rl[,i]))
  print(sum(post_rl[,i] > 0) / length(post_rl[,i]))
}

#shoot length
#set.seed(33)
#model_sl <- brm(Shoot_Length_cm ~ Watering * inoculant, 
#                data = metadata_ep, family = "Gaussian", iter = 2000, 
#                control = list(max_treedepth = 10))
#saveRDS(model_sl, "model_sl.rds")
model_sl <- readRDS("model_sl.rds")
summary(model_sl)
mcmc_plot(model_sl)
plot(conditional_effects(model_sl))
post_sl <- as.data.frame(model_sl) # extract posterior samples from model fit
for(i in 2:ncol(post_sl)){
  print(colnames(post_sl[i]))
  print(sum(post_sl[,i] < 0) / length(post_sl[,i]))
  print(sum(post_sl[,i] > 0) / length(post_sl[,i]))
}
#wet weight
#set.seed(33)
#model_ww <- brm(Whole_Plant_wet_weightg ~ Watering * inoculant, 
#                data = metadata_ep, family = "Gaussian", iter = 2000, 
#                control = list(max_treedepth = 10))
#saveRDS(model_ww, "model_ww.rds")
model_ww <- readRDS("model_ww.rds")
summary(model_ww)
mcmc_plot(model_ww)

#dry weight
boxplot(Dry_Weight ~ Treatment, data = metadata_ep, main = "Whole plant dry weight by Treatment Group")
stripchart(Dry_Weight ~ Treatment, data = metadata_ep, add = TRUE, vertical = TRUE, pch = 1, method = "jitter")
#set.seed(33)
model_wdw <- brm(Dry_Weight ~ Watering * inoculant, 
                 data = metadata_ep, family = "Gaussian", iter = 10000, 
                 control = list(max_treedepth = 10))
#saveRDS(model_wdw, "model_wdw.rds")
model_wdw <- readRDS("model_wdw.rds")
summary(model_wdw)
mcmc_plot(model_wdw)
conditional_effects(model_wdw)
post_wdw <- as.data.frame(model_wdw) # extract posterior samples from model fit
hist(post_wdw$b_inoculantb.amy)
for(i in 2:ncol(post_wdw)){
  print(colnames(post_wdw[i]))
  print(sum(post_wdw[,i] < 0) / length(post_wdw[,i]))
  print(sum(post_wdw[,i] > 0) / length(post_wdw[,i]))
}


#root wet weight
#set.seed(33)
#model_rww <- brm(root_weight_wetg ~ Watering * inoculant, 
#                 data = metadata_ep, family = "Gaussian", iter = 10000, 
#                 control = list(max_treedepth = 10))
#saveRDS(model_rww, "model_rww.rds")
model_rww <- readRDS("model_rww.rds")
summary(model_rww)
mcmc_plot(model_rww)

#root dry weight
#set.seed(33)
#model_rdw <- brm(root_dry_g ~ Watering * inoculant, 
#                 data = metadata_ep, family = "Gaussian", iter = 10000, 
#                 control = list(max_treedepth = 10))
#saveRDS(model_rdw, "model_rdw.rds")
model_rdw <- readRDS("model_rdw.rds")
summary(model_rdw)
mcmc_plot(model_rdw)

#stem wet weight
set.seed(33)
#model_sww <- brm(stem_wetg ~ Watering * inoculant, 
#                 data = metadata_ep, family = "Gaussian", iter = 10000, 
#                 control = list(max_treedepth = 10))
#saveRDS(model_sww, "model_sww.rds")
model_sww <- readRDS("model_sww.rds")
summary(model_sww)
mcmc_plot(model_sww)

#stem dry weight
#set.seed(33)
#model_sdw <- brm(stem_dryg ~ Watering * inoculant, 
#                 data = metadata_ep, family = "Gaussian", iter = 10000, 
#                 control = list(max_treedepth = 10))
#saveRDS(model_sdw, "model_sdw.rds")
model_sdw <- readRDS("model_sdw.rds")
summary(model_sdw)
mcmc_plot(model_sdw)

#leaves wet weight
set.seed(33)
#model_lww <- brm(leaves_wet_g ~ Watering * inoculant, 
#                 data = metadata_ep, family = "Gaussian", iter = 10000, 
#                 control = list(max_treedepth = 10))
#saveRDS(model_lww, "model_lww.rds")
model_lww <- readRDS("model_lww.rds")
summary(model_lww)
mcmc_plot(model_lww)

#leaves dry weight
#set.seed(33)
#model_ldw <- brm(leaves_dry_g ~ Watering * inoculant, 
#                 data = metadata_ep, family = "Gaussian", iter = 10000, 
#                 control = list(max_treedepth = 10))
#saveRDS(model_ldw, "model_ldw.rds")
model_ldw <- readRDS("model_ldw.rds")
summary(model_ldw)
mcmc_plot(model_ldw, main = "effect size distributions for dry weight of leaves")

#%C
boxplot(Weight_percC ~ Treatment, data = metadata_ep, main = "%C by Treatment Group")
stripchart(Weight_percC ~ Treatment, data = metadata_ep, add = TRUE, vertical = TRUE, pch = 1, method = "jitter")
#set.seed(33)
#model_percC <- brm(Weight_percC*.01 ~ Watering * inoculant, 
#                 data = metadata_ep, family = "beta", iter = 10000, 
#                 control = list(max_treedepth = 10))
#saveRDS(model_percC, "model_percC.rds")
model_percC <- readRDS("model_percC.rds")
summary(model_percC)
mcmc_plot(model_percC)
plot(conditional_effects(model_percC))
post_percC <- as.data.frame(model_percC) # extract posterior samples from model fit
hist(post_percC$b_Wateringdrought)
for(i in 2:ncol(post_percC)){
  print(colnames(post_percC[i]))
  print(sum(post_percC[,i] < 0) / length(post_percC[,i]))
  print(sum(post_percC[,i] > 0) / length(post_percC[,i]))
}

#13C
boxplot(Delta13C ~ Treatment, data = metadata_ep)
stripchart(Delta13C ~ Treatment, data = metadata_ep, add = TRUE, vertical = TRUE, pch = 1, method = "jitter")
#set.seed(33)
#model_13C <- brm(Delta13C ~ Watering * inoculant, 
#                 data = metadata_ep, family = "Gaussian", iter = 10000, 
#                 control = list(max_treedepth = 10))
#saveRDS(model_13C, "model_13C.rds")
model_13C <- readRDS("model_13C.rds")
summary(model_13C)
mcmc_plot(model_13C)
plot(conditional_effects(model_13C))

post_13C <- as.data.frame(model_13C) # extract posterior samples from model fit
hist(post_13C$b_Wateringdrought)
for(i in 2:ncol(post_13C)){
  print(colnames(post_13C[i]))
  print(sum(post_13C[,i] < 0) / length(post_13C[,i]))
  print(sum(post_13C[,i] > 0) / length(post_13C[,i]))
}


#C:N
set.seed(33)
model_CN <- brm(CNratio ~ Watering * inoculant, 
                 data = metadata_ep, family = "gamma", iter = 10000, 
                 control = list(max_treedepth = 10))
#saveRDS(model_CN, "model_CN.rds")
summary(model_CN)
plot(conditional_effects(model_CN))


#%N
boxplot(Weight_percN ~ Treatment, data = metadata_ep)
stripchart(Weight_percN ~ Treatment, data = metadata_ep, add = TRUE, vertical = TRUE, pch = 1, method = "jitter")
#set.seed(33)
#model_percN <- brm(Weight_percN*.01 ~ Watering * inoculant, 
#                 data = metadata_ep, family = "beta", iter = 10000, 
#                 control = list(max_treedepth = 10))
#saveRDS(model_percN, "model_percN.rds")
model_percN <- readRDS("model_percN.rds")
summary(model_percN)
mcmc_plot(model_percN)
plot(conditional_effects(model_percN), ylab = "Weight %N")
post_percN <- as.data.frame(model_percN) # extract posterior samples from model fit
hist(post_percN$b_Wateringdrought)
for(i in 2:ncol(post_percN)){
  print(colnames(post_percC[i]))
  print(sum(post_percN[,i] < 0) / length(post_percN[,i]))
  print(sum(post_percN[,i] > 0) / length(post_percN[,i]))
}

#15N
#set.seed(33)
#model_15N <- brm(Delta15N ~ Watering * inoculant, 
#                 data = metadata_ep, family = "Gaussian", iter = 10000, 
#                 control = list(max_treedepth = 10))
#saveRDS(model_15N, "model_15N.rds")
model_15N <- readRDS("model_15N.rds")
summary(model_15N)
mcmc_plot(model_15N)
plot(conditional_effects(model_15N))
post_15N <- as.data.frame(model_15N) # extract posterior samples from model fit
for(i in 2:ncol(post_15N)){
  print(colnames(post_15N[i]))
  print(sum(post_15N[,i] < 0) / length(post_15N[,i]))
  print(sum(post_15N[,i] > 0) / length(post_15N[,i]))
}

#16s Richness
set.seed(33)
boxplot(rich_16s ~ Treatment, data = metadata_ep, main = "16s Richness by Treatment Group")
stripchart(rich_16s ~ Treatment, data = metadata_ep, add = TRUE, vertical = TRUE, pch = 1, method = "jitter")
model_16sR <- brm(rich_16s ~ Watering * inoculant, 
                 data = metadata_ep, family = "negbinomial", iter = 10000, 
                 control = list(max_treedepth = 10))
par(mar = c(5, 4, 4, 2))

boxplot(abund_cult ~ Treatment, data = metadata, main = "Culture abundance by Treatment Group")
stripchart(abund_cult ~ Treatment, data = metadata, add = TRUE, vertical = TRUE, pch = 1, method = "jitter")

model_acult <- brm(abund_cult ~ Watering * inoculant, 
                  data = metadata, family = "negbinomial", iter = 10000, 
                  control = list(max_treedepth = 10))
summary(model_acult)
mcmc_plot(model_acult)

saveRDS(model_16sR, "model_16sR.rds")
model_16sR <- readRDS("model_16sR.rds")
summary(model_16sR)
mcmc_plot(model_16sR)
plot(conditional_effects(model_16sR))

#16s ENS
set.seed(33)
boxplot(ENS_16s ~ Treatment, data = metadata_ep, main = "16s Richness by Treatment Group")
stripchart(ENS_16s ~ Treatment, data = metadata_ep, add = TRUE, vertical = TRUE, pch = 1, method = "jitter")
model_16sENS <- brm(ENS_16s ~ Watering * inoculant, 
                  data = metadata_ep, family = "negbinomial", iter = 10000, 
                  control = list(max_treedepth = 10))
saveRDS(model_16sENS, "model_16sENS.rds")
model_16sENS <- readRDS("model_16sENS.rds")
summary(model_16sENS)
mcmc_plot(model_16sENS)
plot(conditional_effects(model_16sENS))

#ITS Richness
set.seed(33)
boxplot(rich_ITS ~ Treatment, data = metadata_ep, main = "ITS Richness by Treatment Group")
stripchart(rich_ITS ~ Treatment, data = metadata_ep, add = TRUE, vertical = TRUE, pch = 1, method = "jitter")
model_ITSR <- brm(rich_ITS ~ Watering * inoculant, 
                  data = metadata_ep, family = "negbinomial", iter = 10000, 
                  control = list(max_treedepth = 10))
saveRDS(model_ITSR, "model_ITSR.rds")
model_ITSR <- readRDS("model_ITSR.rds")
summary(model_ITSR)
mcmc_plot(model_ITSR)
plot(conditional_effects(model_ITSR))

#ITS ENS
set.seed(33)
boxplot(ENS_ITS ~ Treatment, data = metadata_ep, main = "ITS Richness by Treatment Group")
stripchart(ENS_ITS ~ Treatment, data = metadata_ep, add = TRUE, vertical = TRUE, pch = 1, method = "jitter")
model_ITSENS <- brm(ENS_ITS ~ Watering * inoculant, 
                    data = metadata_ep, family = "negbinomial", iter = 10000, 
                    control = list(max_treedepth = 10))
saveRDS(model_ITSENS, "model_ITSENS.rds")
model_ITSENS <- readRDS("model_ITSENS.rds")
summary(model_ITSENS)
mcmc_plot(model_ITSENS)
plot(conditional_effects(model_ITSENS))

#Beta diversity model
#Calculate bray-curtis
dist.16s <- vegdist(data_16s, method = "bray", na.rm = TRUE)
#Run dbRDA beta diversity model
dbrda_16s <- dbrda(dist.16s ~ Watering * inoculant, data = metadata_ep, na.action = na.exclude, dist = "bray")
summary(dbrda_16s)
#Create vector specifying pch by treatment group (for figure)
pch_type <- as.numeric(lapply(substr(row.names(sites_16s),1,2), function(x) {
  a <- gsub("DB", "1", x)
  b <- gsub("DN", "2", a)
  c <- gsub("DS", "3", b)
  d <- gsub("DW", "4", c)
  e <- gsub("WB", "5", d)
  f <- gsub("WN", "6", e)
  g <- gsub("WS", "7", f)
  h <- gsub("WW", "8", g)
}))
#extract important dbRDA things for figures
sum_16s <- summary(dbrda_16s) #my dbrda is stored asc dbrda_16s
sites_16s <- as.data.frame(sum_16s$sites) #extract oordinates for samples
par(mar = c(5, 4, 4, 2)) #setting figure par
plot(dbRDA2 ~ dbRDA1, data = sites_16s, pch = pch_type, type = "p", cex = 1.5, col = "black", main = "16S beta diversity") #pch_type is a vector of numbers corresponding with my treatment groups
arrows(x0 = rep(0,7), y0 = rep(0,7), x1 = dbrda_16s$CCA$biplot[,1]*1.5, y1 = dbrda_16s$CCA$biplot[,2]*1.5, lwd = 3, col = adjustcolor("black", alpha = .5)) #arrows draws from x0 and y0 to x1 and y1, can make as many arrows as you want as long as you supply a x0,y0,x1,y1 for each arrow.
text(x = dbrda_16s$CCA$biplot[,1]*1.5, y = dbrda_16s$CCA$biplot[,2]*1.5, labels = row.names(dbrda_16s$CCA$biplot)) #text works similarly to arrows, as many different text elements as you want as long as you supply x and y for each
legend(x = "bottomleft", legend= c("DB","DN","DS","DW","WB","WN","WS","WW"), pch = 1:8, pt.cex = 1.5)

#Run permutational anova on dbRDA model
permanova_16S <- anova.cca(dbrda_16s, by = "terms", permutations = 999)
permanova_16S

#nmds
nmds_16s <- metaMDS(dist.16s distance = "bray", iter = 10000, k = 2)
plot(nmds_16s$points, pch = c(1:8)[as.factor(metadata_ep$Treatment)], col = rainbow(8)[as.factor(metadata_ep$Treatment)], main = "16S NMDS")
ordispider(nmds_16s, display = "sites", groups = metadata_ep$Treatment, col = "grey70")
legend(x = "bottomright", legend= c("DB","DN","DS","DW","WB","WN","WS","WW"), pch = 1:8, pt.cex = 1.5)

#Calculate bray-curtis
dist.ITS <- vegdist(data_ITS, method = "bray", na.rm = TRUE)
#Run dbRDA beta diversity model
dbrda_ITS <- dbrda(dist.ITS ~ Watering * inoculant, data = metadata_ep, na.action = na.exclude, dist = "bray")

#extract important dbRDA things for figures
sum_ITS <- summary(dbrda_ITS)
biplot <- sum_ITS_biplot
sites_ITS <- as.data.frame(sum_ITS$sites)
plot(dbRDA2 ~ dbRDA1, data = sites_ITS, pch = pch_type, type = "p", cex = 1.5, col = "black", main = "ITS beta diversity")
arrows(x0 = rep(0,7), y0 = rep(0,7), x1 = dbrda_ITS$CCA$biplot[,1]*1.5, y1 = dbrda_ITS$CCA$biplot[,2]*1.5, lwd = 3, col = adjustcolor("black", alpha = .5))
text(x = dbrda_ITS$CCA$biplot[,1]*1.5, y = dbrda_ITS$CCA$biplot[,2]*1.5, labels = row.names(dbrda_ITS$CCA$biplot))
legend(x = "topright", legend= c("DB","DN","DS","DW","WB","WN","WS","WW"), pch = 1:8, pt.cex = 1.5)

#Run permutational anova on dbRDA model
permanova_ITS <- anova.cca(dbrda_ITS, by = "terms", permutations = 999)
permanova_ITS

#nmds
nmds_ITS <- metaMDS(dist.ITS, distance = "bray", iter = 10000, k = 2)
plot(nmds_ITS$points, pch = c(1:8)[as.factor(metadata_ep$Treatment)], col = rainbow(8)[as.factor(metadata_ep$Treatment)], main = "ITS NMDS")
ordispider(nmds_ITS, display = "sites", groups = metadata_ep$Treatment, col = "grey70")
legend(x = "topright", legend= c("DB","DN","DS","DW","WB","WN","WS","WW"), pch = 1:8, pt.cex = 1.5)

#Misc. Figures ####
plot(Photo ~ var_photo, data = metadata[c(1:49,51:480),])
abline(lm(Photo ~ var_photo, data = metadata[c(1:49,51:480),]), col = "red")
summary(lm(Photo ~ var_photo, data = metadata[c(1:49,51:480),]))
boxplot(var_photo ~ Treatment, data = metadata[c(1:49,51:480),])
boxplot(Photo ~ Treatment, data = metadata[c(1:49,51:480),])
glm(Treament)
#Wbj1, row 50, extreme photo outlier
palette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
#Photosynthesis
plot(Photo ~ Date_cont, data = metadata[c(1:49,51:480),], 
     col = adjustcolor(palette, alpha = .75)[as.factor(metadata$Treatment)], pch = 16,
     ylab = "photosynthesis", xlab = "time point", main = "Photosynthesis at each time point")
abline(lm(Photo ~ Date_cont, data = metadata[c(1:49,51:480),][metadata$Treatment == "DB",]), 
       col = adjustcolor(palette[1], alpha = .5), lwd = 6)
text(x = 6, y = -1, labels = "DB", col = palette[1], cex = 1.5)
abline(lm(Photo ~ Date_cont, data = metadata[c(1:49,51:480),][metadata$Treatment == "DN",]), 
       col = adjustcolor(palette[2], alpha = .5), lwd = 6)
text(x = 6, y = 15.5, labels = "DN", col = palette[2], cex = 1.5)
abline(lm(Photo ~ Date_cont, data = metadata[c(1:49,51:480),][metadata$Treatment == "DS",]), 
       col = adjustcolor(palette[3], alpha = .5), lwd = 6)
text(x = 6, y = 7, labels = "DS", col = palette[3], cex = 1.5)
abline(lm(Photo ~ Date_cont, data = metadata[c(1:49,51:480),][metadata$Treatment == "DW",]),
       col = adjustcolor(palette[4], alpha = .5), lwd = 6)
text(x = 5.5, y = 4, labels = "DW", col = palette[4], cex = 1.5)
abline(lm(Photo ~ Date_cont, data = metadata[c(1:49,51:480),][metadata$Treatment == "WB",]), 
       col = adjustcolor(palette[5], alpha = .5), lwd = 6)
text(x = 6, y = 3.5, labels = "WB", col = palette[5], cex = 1.5)
abline(lm(Photo ~ Date_cont, data = metadata[c(1:49,51:480),][metadata$Treatment == "WN",]), 
       col = adjustcolor(palette[6], alpha = .5), lwd = 6)
text(x = 6, y = 26, labels = "WN", col = palette[6], cex = 1.5)
abline(lm(Photo ~ Date_cont, data = metadata[c(1:49,51:480),][metadata$Treatment == "WS",]), 
       col = adjustcolor(palette[7], alpha = .5), lwd = 6)
text(x = 6, y = 20.5, labels = "WS", col = palette[7], cex = 1.5)
abline(lm(Photo ~ Date_cont, data = metadata[c(1:49,51:480),][metadata$Treatment == "WW",]), 
       col = adjustcolor(palette[8], alpha = .5), lwd = 6)
text(x = 5.5, y = -1, labels = "WW", col = palette[8], cex = 1.5)

#Variance from photosynthesis mean
plot(var_photo ~ Date_cont, data = metadata[c(1:49,51:480),], 
     col = adjustcolor(palette, alpha = .75)[as.factor(metadata$Treatment)], pch = 16,
     ylab = "photo - mean(photo)", xlab = "time point", main = "Variance from mean photosynthesis at each time point")
abline(lm(var_photo ~ Date_cont, data = metadata[c(1:49,51:480),][metadata$Treatment == "DB",]), 
       col = adjustcolor(palette[1], alpha = .5), lwd = 6)
text(x = 6, y = -5, labels = "DB", col = palette[1], cex = 1.5)
abline(lm(var_photo ~ Date_cont, data = metadata[c(1:49,51:480),][metadata$Treatment == "DN",]), 
       col = adjustcolor(palette[2], alpha = .5), lwd = 6)
text(x = 6, y = 10.5, labels = "DN", col = palette[2], cex = 1.5)
abline(lm(var_photo ~ Date_cont, data = metadata[c(1:49,51:480),][metadata$Treatment == "DS",]), 
       col = adjustcolor(palette[3], alpha = .5), lwd = 6)
text(x = 5.5, y = 0, labels = "DS", col = palette[3], cex = 1.5)
abline(lm(var_photo ~ Date_cont, data = metadata[c(1:49,51:480),][metadata$Treatment == "DW",]),
       col = adjustcolor(palette[4], alpha = .5), lwd = 6)
text(x = 5.5, y = -4, labels = "DW", col = palette[4], cex = 1.5)
abline(lm(var_photo ~ Date_cont, data = metadata[c(1:49,51:480),][metadata$Treatment == "WB",]), 
       col = adjustcolor(palette[5], alpha = .5), lwd = 6)
text(x = 6, y = -1, labels = "WB", col = palette[5], cex = 1.5)
abline(lm(var_photo ~ Date_cont, data = metadata[c(1:49,51:480),][metadata$Treatment == "WN",]), 
       col = adjustcolor(palette[6], alpha = .5), lwd = 6)
text(x = 6, y = 13, labels = "WN", col = palette[6], cex = 1.5)
abline(lm(var_photo ~ Date_cont, data = metadata[c(1:49,51:480),][metadata$Treatment == "WS",]), 
       col = adjustcolor(palette[7], alpha = .5), lwd = 6)
text(x = 6, y = 16, labels = "WS", col = palette[7], cex = 1.5)
abline(lm(var_photo ~ Date_cont, data = metadata[c(1:49,51:480),][metadata$Treatment == "WW",]), 
       col = adjustcolor(palette[8], alpha = .5), lwd = 6)
text(x = 6, y = -9, labels = "WW", col = palette[8], cex = 1.5)

#Canopy area
plot(`Canopy Area (cm^2)` ~ Date_cont, data = metadata[c(1:49,51:480),], 
     col = adjustcolor(palette, alpha = .75)[as.factor(metadata$Treatment)], pch = 16,
     ylab = "Canopy Area (cm^2)", xlab = "time point", main = "Canopy Area at each time point")
abline(lm(`Canopy Area (cm^2)` ~ Date_cont, data = metadata[c(1:49,51:480),][metadata$Treatment == "DB",]), 
       col = adjustcolor(palette[1], alpha = .5), lwd = 6)
text(x = 6, y = -5, labels = "DB", col = palette[1], cex = 1.5)
abline(lm(`Canopy Area (cm^2)` ~ Date_cont, data = metadata[c(1:49,51:480),][metadata$Treatment == "DN",]), 
       col = adjustcolor(palette[2], alpha = .5), lwd = 6)
text(x = 6, y = 10.5, labels = "DN", col = palette[2], cex = 1.5)
abline(lm(`Canopy Area (cm^2)` ~ Date_cont, data = metadata[c(1:49,51:480),][metadata$Treatment == "DS",]), 
       col = adjustcolor(palette[3], alpha = .5), lwd = 6)
text(x = 5.5, y = 0, labels = "DS", col = palette[3], cex = 1.5)
abline(lm(`Canopy Area (cm^2)` ~ Date_cont, data = metadata[c(1:49,51:480),][metadata$Treatment == "DW",]),
       col = adjustcolor(palette[4], alpha = .5), lwd = 6)
text(x = 5.5, y = -4, labels = "DW", col = palette[4], cex = 1.5)
abline(lm(`Canopy Area (cm^2)` ~ Date_cont, data = metadata[c(1:49,51:480),][metadata$Treatment == "WB",]), 
       col = adjustcolor(palette[5], alpha = .5), lwd = 6)
text(x = 6, y = -1, labels = "WB", col = palette[5], cex = 1.5)
abline(lm(`Canopy Area (cm^2)` ~ Date_cont, data = metadata[c(1:49,51:480),][metadata$Treatment == "WN",]), 
       col = adjustcolor(palette[6], alpha = .5), lwd = 6)
text(x = 6, y = 13, labels = "WN", col = palette[6], cex = 1.5)
abline(lm(`Canopy Area (cm^2)` ~ Date_cont, data = metadata[c(1:49,51:480),][metadata$Treatment == "WS",]), 
       col = adjustcolor(palette[7], alpha = .5), lwd = 6)
text(x = 6, y = 16, labels = "WS", col = palette[7], cex = 1.5)
abline(lm(`Canopy Area (cm^2)` ~ Date_cont, data = metadata[c(1:49,51:480),][metadata$Treatment == "WW",]), 
       col = adjustcolor(palette[8], alpha = .5), lwd = 6)
text(x = 6, y = -9, labels = "WW", col = palette[8], cex = 1.5)

#canopy area
boxplot(`Canopy_Area_cm^2` ~ Date, data = metadata[metadata$Group == "WW",], na.rm = TRUE, col = "Blue", ylim = c(0,20))
boxplot(`Canopy_Area_cm^2` ~ Date, data = metadata[metadata$Group == "WS",], na.rm = TRUE, col = "Red", add = TRUE)

CA_mean <- tapply(X = metadata$`Canopy_Area_cm^2`, INDEX = list(metadata$Date, metadata$Group), FUN = mean, na.rm = TRUE)
CA_sd <- tapply(X = metadata$`Canopy_Area_cm^2`, INDEX = list(metadata$Date, metadata$Group), FUN = sd, na.rm = TRUE)
CA_graph <- melt(CA_mean)
colnames(CA_graph) <- c("Date", "Group", "mean_CA")
CA_graph <- cbind(CA_graph, "d_lin" = rep(c(1:6), 8), "SD" = melt(CA_sd)$value)
CA_graph <- cbind(CA_graph, "error_plus" = c(CA_graph$mean_CA + CA_graph$SD), "error_minus" = CA_graph$mean_CA - CA_graph$SD)

par(mar = c(5, 4, 4, 2), bty = "n")
plot(mean_CA ~ d_lin, data = CA_graph[CA_graph$Group == "WW",], type = "l", lwd = 7, col = alpha("seagreen4", .25), ylim = c(0,17), ylab = "Average Canopy Area (cm^2)", xlab = "Week", main = "Change in canopy area over time")
points(mean_CA ~ d_lin, data = CA_graph[CA_graph$Group == "WW",], type = "p", col = alpha("seagreen4", 1), pch = 19)
#mtext(text = "whole community, watered", side = 4, las = 1, at = 4.278700, col = "seagreen4")
text(x = 6, y = 4, labels = "whole community, watered", pos = 2, col = "seagreen4")
text(x = 6, y = 3.5, labels = "*p < 0.01", pos = 2, col = "seagreen4")
#arrows(x0 = c(1:6), y0 = CA_graph[CA_graph$Group == "WW",]$error_plus, 
#       x1 = c(1:6), y1 = CA_graph[CA_graph$Group == "WW",]$error_minus,
#       angle = 90, col = alpha("red", .25), lwd = 2, length = .075, code = 3)
points(mean_CA ~ d_lin, data = CA_graph[CA_graph$Group == "WS",], type = "l", lwd = 7, col = alpha("darkblue", .25))
points(mean_CA ~ d_lin, data = CA_graph[CA_graph$Group == "WS",], type = "p", col = alpha("darkblue", 1), pch = 19)
text(x = 6, y = 13.5, labels = "sterile water, watered", pos = 2, col = "darkblue")
text(x = 6, y = 13, labels = "*p = 0.01", pos = 2, col = "darkblue")
#arrows(x0 = c(1:6), y0 = CA_graph[CA_graph$Group == "WS",]$error_plus, 
#       x1 = c(1:6), y1 = CA_graph[CA_graph$Group == "WS",]$error_minus,
#       angle = 90, col = alpha("blue", .25), lwd = 2, length = .075, code = 3)
points(mean_CA ~ d_lin, data = CA_graph[CA_graph$Group == "DS",], type = "l", lwd = 5, col = alpha("royalblue", .25))
points(mean_CA ~ d_lin, data = CA_graph[CA_graph$Group == "DS",], type = "p", col = alpha("royalblue", .5), pch = 19)
text(x = 6, y = 10.5, labels = "sterile water, drought", pos = 2, col = "royalblue")
text(x = 6, y = 10, labels = "p = 0.39", pos = 2, col = "royalblue")
#arrows(x0 = c(1:6), y0 = CA_graph[CA_graph$Group == "DS",]$error_plus, 
#       x1 = c(1:6), y1 = CA_graph[CA_graph$Group == "DS",]$error_minus,
#       angle = 90, col = alpha("blue", .25), lwd = 2, length = .075, code = 3)
points(mean_CA ~ d_lin, data = CA_graph[CA_graph$Group == "DW",], type = "l", lwd = 5, col = alpha("seagreen2", .25))
points(mean_CA ~ d_lin, data = CA_graph[CA_graph$Group == "DW",], type = "p", col = alpha("seagreen2", .5), pch = 19)
text(x = 6, y = 8, labels = "whole community, drought", pos = 2, col = "seagreen2")
text(x = 6, y = 7.5, labels = "p = 0.28", pos = 2, col = "seagreen2")
#arrows(x0 = c(1:6), y0 = CA_graph[CA_graph$Group == "DW",]$error_plus, 
#       x1 = c(1:6), y1 = CA_graph[CA_graph$Group == "DW",]$error_minus,
#       angle = 90, col = alpha("blue", .25), lwd = 2, length = .075, code = 3)

WW_CA_m <- glmer.nb(`Canopy_Area_cm^2` ~ Date_cont + (1|Plant_ID), data = metadata[metadata$Group == "WW",])
summary(WW_CA_m)
WS_CA_m <- glmer.nb(`Canopy_Area_cm^2` ~ Date_cont + (1|Plant_ID), data = metadata[metadata$Group == "WS",])
summary(WS_CA_m)
DW_CA_m <- glmer.nb(`Canopy_Area_cm^2` ~ Date_cont + (1|Plant_ID), data = metadata[metadata$Group == "DW",])
summary(DW_CA_m)
DS_CA_m <- glmer.nb(`Canopy_Area_cm^2` ~ Date_cont + (1|Plant_ID), data = metadata[metadata$Group == "DS",])
summary(DS_CA_m)
WB_CA_m <- glmer.nb(`Canopy_Area_cm^2` ~ Date_cont + (1|Plant_ID), data = metadata[metadata$Group == "WB",])
summary(WB_CA_m)
DB_CA_m <- glmer.nb(`Canopy_Area_cm^2` ~ Date_cont + (1|Plant_ID), data = metadata[metadata$Group == "DB",])
summary(DB_CA_m)
WN_CA_m <- glmer.nb(`Canopy_Area_cm^2` ~ Date_cont + (1|Plant_ID), data = metadata[metadata$Group == "WN",])
summary(WN_CA_m)
DN_CA_m <- glmer.nb(`Canopy_Area_cm^2` ~ Date_cont + (1|Plant_ID), data = metadata[metadata$Group == "DN",])
summary(DN_CA_m)




#DNA
metadata_last <- metadata[metadata$Date == "6/30/2023",]

data_16s <- data_16s[row.names(data_16s) %in% unique(metadata_last$Plant_ID),]
data_ITS <- data_ITS[row.names(data_ITS) %in% unique(metadata_last$Plant_ID),]

colnames(metadata_last)
metadata_last <- metadata_last[,1:25]

summary(rowSums(data_16s))
sort(rowSums(data_16s))
summary(rowSums(data_ITS))

set.seed(123)
data_16s.r <- rrarefy(x = data_16s, sample = 1000)
set.seed(123)
data_ITS.r <- rrarefy(x = data_ITS, sample = 2910)

metadata_last <- cbind(metadata_last, "Richness_16s" = rowSums(data_16s.r > 0))
metadata_last <- cbind(metadata_last, "Shannon_16s" = diversity(data_16s.r))
metadata_last <- cbind(metadata_last, "Richness_ITS" = rowSums(data_ITS.r > 0))
metadata_last <- cbind(metadata_last, "Shannon_ITS" = diversity(data_ITS.r))

hist(rowSums(data_16s.r > 0))
hist(rowSums(data_ITS.r > 0))

boxplot(metadata_last$Shannon_16s ~ metadata_last$inoculant)
t.test(metadata_last$Richness_16s ~ metadata_last$inoculant)

boxplot(metadata_last$Shannon_ITS ~ metadata_last$inoculant)
t.test(metadata_last$Richness_ITS ~ metadata_last$Watering)



plot(metadata_last$Richness_16s ~ metadata_last$`Dry Weight`, ylab = "16s Richness", xlab = "Pot Weight (g)", main = "Richness across pot weights")
abline(lm(metadata_last$Richness_16s ~ metadata_last$`Dry Weight`), col = "red")
cor.test(metadata_last$Richness_16s, metadata_last$`Dry Weight`)

plot(metadata_last$Richness_ITS ~ metadata_last$`Dry Weight`, ylab = "ITS Richness", xlab = "Pot Weight (g)", main = "Richness across pot werights")
abline(lm(metadata_last$Richness_ITS ~ metadata_last$`Dry Weight`), col = "red")
cor.test(metadata_last$Richness_ITS, metadata_last$`Dry Weight`)

plot(metadata_last$Shannon_16s ~ metadata_last$`Pot Weight`, ylab = "16s Shannon", xlab = "Pot Weight (g)", main = "Shannon across pot werights")
abline(lm(metadata_last$Shannon_16s ~ metadata_last$`Pot Weight`), col = "red")
cor.test(metadata_last$Shannon_16s, metadata_last$`Pot Weight`)

plot(metadata_last$Shannon_ITS ~ metadata_last$`Pot Weight`, ylab = "ITS Shannon", xlab = "Pot Weight (g)", main = "Shannon across pot werights")
abline(lm(metadata_last$Shannon_ITS ~ metadata_last$`Pot Weight`), col = "red")
cor.test(metadata_last$Shannon_ITS, metadata_last$`Pot Weight`)

m1 <- brm(Shannon_16s ~ Watering * inoculant + (1|Plant_ID), data = metadata_last, iter = 20000, control = list(adapt_delta = 0.9, max_treedepth = 20), family = Gamma(link = "log"))

summary(m1)
mcmc_plot(m1)

m2 <- brm(Shannon_ITS ~ Watering * inoculant + (1|Plant_ID), data = metadata_last, iter = 20000, control = list(adapt_delta = 0.9, max_treedepth = 20), family = Gamma(link = "log"))

summary(m2)

pc_data <- data.frame(metadata_last[,c(6:16,24:25,26:29)])
prcomp(as.numeric(pc_data))

is.infinite(pc_data)


nmds_16s <- metaMDS(data_16s.r, try = 200)
ordiplot(nmds_16s, type = "t", display = "sites")
ordihull(nmds_16s, metadata_last$Watering)

nmds_ITS <- metaMDS(data_ITS.r, try = 200)
ordiplot(nmds_ITS, type = "t", display = "sites")
ordihull(nmds_ITS, metadata_last$inoculant)

# subset taxonomy by species
#Collapse data by taxonomy
genus_ITS <- data.frame("genus" = f_tax$genus, t(data_ITS))
genus_ITS <- aggregate(. ~ genus_ITS[,1], genus_ITS[,2:ncol(genus_ITS)], sum)

data_prop <- c()
#genus_ITS <- genus_ITS[-c(1,10,26),]
for (x in 2:ncol(genus_ITS)){
  data_prop <- cbind(data_prop, genus_ITS[,x] / sum(genus_ITS[,x]))}

row.names(data_prop) <- genus_ITS[,1]
colnames(data_prop) <- colnames(genus_ITS[,2:ncol(genus_ITS)])
data_prop <- data_prop[order(rowSums(data_prop[,2:ncol(data_prop)]),decreasing = T),]
data_prop <- rbind(data_prop[1:10,], "other" = colSums(data_prop[11:nrow(data_prop),]))
customcol <- c("cadetblue4","royalblue3","darkblue","tomato1","dodgerblue2",
               "cyan","darkred","purple","mediumblue","palegoldenrod",
               "grey", "lightgoldenrod","indianred","yellow","purple4","darkgreen",
               "lightsalmon","yellow3","purple2","lightblue","firebrick",
               "navy","red4","red","darkmagenta","mediumvioletred",
               "violetred2","skyblue","dodgerblue4")
par(mar = c(5, 1, 4, 11))
barplot(data_prop, col=customcol, main = paste("ITS Taxa barplot at genus level", sep = " "), legend.text = F, axes = F, cex.names = .8, las = 2, border=NA, space=-0.1)
legend(x = "topright", inset = c(-0.2, 0), legend = row.names(data_prop),
       fill = customcol, xpd = TRUE, cex = .75, border = NA, box.col = "white")

# subset taxonomy by species
#Collapse data by taxonomy
genus_16s <- data.frame("genus" = b_tax$genus, t(data_16s))
genus_16s <- aggregate(. ~ genus_16s[,1], genus_16s[,2:ncol(genus_16s)], sum)

data_prop <- c()
for (x in 2:ncol(genus_16s)){
  data_prop <- cbind(data_prop, genus_16s[,x] / sum(genus_16s[,x]))}
row.names(data_prop) <- genus_16s[,1]
colnames(data_prop) <- colnames(genus_16s[,2:ncol(genus_16s)])
data_prop <- data_prop[order(rowSums(data_prop, na.rm = TRUE), decreasing = T),]
data_prop <- rbind(data_prop[1:10,], "other" = colSums(data_prop[11:nrow(data_prop),]))
customcol <- c("cadetblue4","royalblue3","darkblue","tomato1","dodgerblue2",
                "cyan","darkred","purple","mediumblue","palegoldenrod",
               "grey","lightgoldenrod","indianred","yellow","purple4","darkgreen",
               "lightsalmon","yellow3","purple2","lightblue","firebrick",
               "navy","red4","red","darkmagenta","mediumvioletred",
               "violetred2","skyblue","dodgerblue4")
par(mar = c(5, 1, 4, 11))
barplot(data_prop, col=customcol, main = paste("16S Taxa barplot at genus level", sep = " "), legend.text = F, axes = F, cex.names = .8, las = 2, border=NA, space=-0.1)
legend(x = "topright", inset = c(-0.2, 0), legend = row.names(data_prop),
       fill = customcol, xpd = TRUE, cex = .75, border = NA, box.col = "white")

#ANCOM-BC at genus level
#ITS
ITS_taxtable <- cbind(f_tax[,1:6], "species" = paste(row.names(f_tax), "_", f_tax$species, sep = ""))
ITS_taxtable <- tax_table(ITS_taxtable)
row.names(ITS_taxtable) <- row.names(f_tax)
colnames(ITS_taxtable) <- colnames(f_tax)
ITS_phyloseq <- phyloseq(otu_table = otu_table(t(data_ITS), taxa_are_rows = TRUE), tax_table = ITS_taxtable, sample_data(metadata))

results_ITS_genus <- ancombc2(data = ITS_phyloseq, tax_level = "genus", fix_formula = "Watering + inoculant", p_adj_method = "holm", prv_cut = 0.10, group = NULL, struc_zero = FALSE, neg_lb = FALSE, alpha = 0.05)

View(results_ITS_genus$res)

results_ITS_genus$res[results_ITS_genus$res$diff_Wateringdrought == TRUE,]$taxon

##ITS lfc graphs
par(mar = c(5, 4, 4, 2))
#drought
plot(results_ITS_genus$res[results_ITS_genus$res$diff_Wateringdrought == TRUE,]$lfc_Wateringdrought[c(1,2,4)], ylab = "log fold change", main = "ITS drought lfc", ylim = c(-1.6,0), xlim = c(0,4))
text(x = c(1,2,3), y = results_ITS_genus$res[results_ITS_genus$res$diff_Wateringdrought == TRUE,]$lfc_Wateringdrought[c(1,2,4)] + .1, labels = c("Cladosporium", "Martensia", "Cystobasidium"))
arrows(y0 = results_ITS_genus$res[results_ITS_genus$res$diff_Wateringdrought == TRUE,]$lfc_Wateringdrought + results_ITS_genus$res[results_ITS_genus$res$diff_Wateringdrought == TRUE,]$se_Wateringdrought, y1 = results_ITS_genus$res[results_ITS_genus$res$diff_Wateringdrought == TRUE,]$lfc_Wateringdrought - results_ITS_genus$res[results_ITS_genus$res$diff_Wateringdrought == TRUE,]$se_Wateringdrought, x0 = 1:5, x1 = 1:5, angle = 90, code = 3)
#sterile water
plot(results_ITS_genus$res[results_ITS_genus$res$`diff_inoculantsterile water` == TRUE,]$`lfc_inoculantsterile water`, ylab = "log fold change", main = "ITS sterile water lfc (acremonium)", ylim = c(0,5))
arrows(y0 = results_ITS_genus$res[results_ITS_genus$res$`diff_inoculantsterile water` == TRUE,]$`lfc_inoculantsterile water` + results_ITS_genus$res[results_ITS_genus$res$`diff_inoculantsterile water` == TRUE,]$`se_inoculantsterile water`, y1 = results_ITS_genus$res[results_ITS_genus$res$`diff_inoculantsterile water` == TRUE,]$`lfc_inoculantsterile water` - results_ITS_genus$res[results_ITS_genus$res$`diff_inoculantsterile water` == TRUE,]$`se_inoculantsterile water`, x0 = 1:5, x1 = 1:5, angle = 90, code = 3)
#bamy
plot(results_ITS_genus$res[results_ITS_genus$res$diff_inoculantb.amy == TRUE,]$lfc_inoculantb.amy, ylab = "log fold change", main = "ITS bamy lfc(penecillium)", ylim = c(0,2))
arrows(y0 = results_ITS_genus$res[results_ITS_genus$res$diff_inoculantb.amy == TRUE,]$lfc_inoculantb.amy + results_ITS_genus$res[results_ITS_genus$res$diff_inoculantb.amy == TRUE,]$se_inoculantb.amy, y1 = results_ITS_genus$res[results_ITS_genus$res$diff_inoculantb.amy == TRUE,]$lfc_inoculantb.amy - results_ITS_genus$res[results_ITS_genus$res$diff_inoculantb.amy == TRUE,]$se_inoculantb.amy, x0 = 1:5, x1 = 1:5, angle = 90, code = 3)
#whole comm
plot(results_ITS_genus$res[results_ITS_genus$res$`diff_inoculantwhole community` == TRUE,]$`lfc_inoculantwhole community`, ylab = "log fold change", main = "ITS wholecomm lfc(penecillium)", ylim = c(0,2))
arrows(y0 = results_ITS_genus$res[results_ITS_genus$res$`diff_inoculantwhole community` == TRUE,]$`lfc_inoculantwhole community` + results_ITS_genus$res[results_ITS_genus$res$`diff_inoculantwhole community` == TRUE,]$`se_inoculantwhole community`, y1 = results_ITS_genus$res[results_ITS_genus$res$`diff_inoculantwhole community` == TRUE,]$`lfc_inoculantwhole community` - results_ITS_genus$res[results_ITS_genus$res$`diff_inoculantwhole community` == TRUE,]$`se_inoculantwhole community`, x0 = 1:5, x1 = 1:5, angle = 90, code = 3)

#16S
B_taxtable <- tax_table(b_tax)
row.names(B_taxtable) <- row.names(b_tax)
colnames(B_taxtable) <- colnames(b_tax)
B_phyloseq <- phyloseq(otu_table = otu_table(t(data_16s), taxa_are_rows = TRUE), tax_table = B_taxtable, sample_data(metadata))

results_B_genus <- ancombc2(data = B_phyloseq, tax_level = "genus", fix_formula = "Watering + inoculant", p_adj_method = "holm", prv_cut = 0.10, group = NULL, struc_zero = FALSE, neg_lb = FALSE, alpha = 0.05)


View(results_B_genus$res)

##B lfc graphs
par(mar = c(5, 4, 4, 2))
#drought
plot(results_B_genus$res[results_B_genus$res$diff_Wateringdrought == TRUE,]$lfc_Wateringdrought, ylab = "log fold change", main = "16s drought lfc", ylim = c(-2,4))
arrows(y0 = results_B_genus$res[results_B_genus$res$diff_Wateringdrought == TRUE,]$lfc_Wateringdrought + results_B_genus$res[results_B_genus$res$diff_Wateringdrought == TRUE,]$se_Wateringdrought, y1 = results_B_genus$res[results_B_genus$res$diff_Wateringdrought == TRUE,]$lfc_Wateringdrought - results_B_genus$res[results_B_genus$res$diff_Wateringdrought == TRUE,]$se_Wateringdrought, x0 = 1:19, x1 = 1:19, angle = 90, code = 3)
#sterile water
plot(results_B_genus$res[results_B_genus$res$`diff_inoculantsterile water` == TRUE,]$`lfc_inoculantsterile water`, ylab = "log fold change", main = "16s sterile water (Mycobacteriaceae)", ylim = c(-2.5,0))
arrows(y0 = results_B_genus$res[results_B_genus$res$`diff_inoculantsterile water` == TRUE,]$`lfc_inoculantsterile water` + results_B_genus$res[results_B_genus$res$`diff_inoculantsterile water` == TRUE,]$`se_inoculantsterile water`, y1 = results_B_genus$res[results_B_genus$res$`diff_inoculantsterile water` == TRUE,]$`lfc_inoculantsterile water` - results_B_genus$res[results_B_genus$res$`diff_inoculantsterile water` == TRUE,]$`se_inoculantsterile water`, x0 = 1:5, x1 = 1:5, angle = 90, code = 3)
#bamy
plot(results_B_genus$res[results_B_genus$res$diff_inoculantb.amy == TRUE,]$lfc_inoculantb.amy, ylab = "log fold change", main = "16s bamy lfc(bacillus, Oxalobacteraceae)", ylim = c(-3,3))
arrows(y0 = results_B_genus$res[results_B_genus$res$diff_inoculantb.amy == TRUE,]$lfc_inoculantb.amy + results_B_genus$res[results_B_genus$res$diff_inoculantb.amy == TRUE,]$se_inoculantb.amy, y1 = results_B_genus$res[results_B_genus$res$diff_inoculantb.amy == TRUE,]$lfc_inoculantb.amy - results_B_genus$res[results_B_genus$res$diff_inoculantb.amy == TRUE,]$se_inoculantb.amy, x0 = 1:5, x1 = 1:5, angle = 90, code = 3)
#whole comm
plot(results_B_genus$res[results_B_genus$res$`diff_inoculantwhole community` == TRUE,]$`lfc_inoculantwhole community`, ylab = "log fold change", main = "16s wholecomm lfc(Rhizobiaceae, Weeksellaceae)", ylim = c(-3, 0))
arrows(y0 = results_B_genus$res[results_B_genus$res$`diff_inoculantwhole community` == TRUE,]$`lfc_inoculantwhole community` + results_B_genus$res[results_B_genus$res$`diff_inoculantwhole community` == TRUE,]$`se_inoculantwhole community`, y1 = results_B_genus$res[results_B_genus$res$`diff_inoculantwhole community` == TRUE,]$`lfc_inoculantwhole community` - results_B_genus$res[results_B_genus$res$`diff_inoculantwhole community` == TRUE,]$`se_inoculantwhole community`, x0 = 1:5, x1 = 1:5, angle = 90, code = 3)
