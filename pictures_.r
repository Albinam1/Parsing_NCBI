require("tidyverse")

effect_type = c('Transcription_factors', 'Gene_expression', 'Protein_structure', 'Alternative_splicing', 'Gene_regulatory_network', 'Epigenetics')

allDf = NULL
for (effect in effect_type) {
  inFile = paste0(effect,"_Result_by_Year.csv")
  inDf = read.csv(inFile)
  colnames(inDf)[2] <- effect
  if(is.null(allDf)){
    allDf = inDf
  }else{
    allDf = allDf %>% full_join(inDf, by = "Year")
  }
}
rm(inDf) 

data <- allDf %>% pivot_longer(!Year, names_to = "group", values_to = "count")
data$value = log(data$count)

empty_bar <- 3
to_add <- data.frame( matrix(NA, empty_bar*nlevels(data$group), ncol(data)) )
colnames(to_add) <- colnames(data)
to_add$group <- rep(levels(data$group), each=empty_bar)
data <- rbind(data, to_add)
data <- data %>% arrange(match(group, effect_type))
data$id <- seq(1, nrow(data))

label_data <- data
number_of_bar <- nrow(label_data)
angle <- 90 - 360 * (label_data$id-0.5) /number_of_bar    
label_data$hjust <- ifelse( angle < -90, 1, 0)
label_data$angle <- ifelse(angle < -90, angle+180, angle)


base_data <- data %>% 
  group_by(group) %>% 
  summarize(start=min(id), end=max(id)) %>% 
  rowwise() %>% 
  mutate(title=mean(c(start, end)))%>% 
  arrange(match(group, effect_type))

base_angle <- label_data %>% group_by(group) %>% 
  summarize(angle = median(angle))
base_data$angle <- ifelse(base_angle$angle > 0, base_angle$angle-90, base_angle$angle-60)

angle = -360 * base_data$title/nrow(data)    
base_data$angle <- ifelse(angle < -90 , angle+180, ifelse(angle < -180 , angle-180, ifelse(angle < -270 , angle+360, angle)))
base_data$angle <- ifelse(angle < -270 , angle+360, ifelse(angle < -180 , angle-180, ifelse(angle < -90 , angle+180, angle)))


grid_data <- base_data
grid_data$end <- grid_data$end[ c( nrow(grid_data), 1:nrow(grid_data)-1)] + 1
grid_data$start <- grid_data$start - 1
grid_data <- grid_data[-1,]

scales = c(2, 4, 6, 8)

p <- ggplot(data, aes(x=as.factor(id), y=value, fill=group)) +    
  
  geom_bar(aes(x=as.factor(id), y=value, fill=group), stat="identity", alpha=0.5) +
  

  geom_segment(data=grid_data, aes(x = end, y = scales[4], xend = start, yend = scales[4]), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = scales[3], xend = start, yend = scales[3]), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = scales[2], xend = start, yend = scales[2]), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = scales[1], xend = start, yend = scales[1]), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  annotate("text", x = rep(max(data$id),4), y = scales, label = scales , color="grey", size=2.5 , angle=0, fontface="bold", hjust=1) +
  annotate("text", x = max(data$id)*0.985, y = (max(scales)+min(scales))/2, label = "log2(count)" , color="grey", size=2.5 , angle=95, fontface="bold", hjust=0.5) +
  
  geom_bar(aes(x=as.factor(id), y=value, fill=group), stat="identity", alpha=0.5) +
  ylim(-min(scales)*4,max(scales)+min(scales)) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-1,4), "cm") 
  )+
  coord_polar() + 
  geom_text(data=label_data, aes(x=id, y=value+min(scales)/5, label=count, hjust=hjust), 
            color="black", alpha=0.6, size=2, angle= label_data$angle, inherit.aes = FALSE ) +
  geom_text(data=label_data, aes(x=id, y=-min(scales)/2, label=Year, hjust=hjust), 
            color="black", alpha=0.6, size=1.5, angle= label_data$angle, inherit.aes = FALSE ) +
  geom_segment(data=base_data, aes(x = start, y = -min(scales)*0.7, xend = end, yend = -min(scales)*0.7), colour = "black", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_text(data=base_data, aes(x = title, y = -min(scales)*1.1 , label=group), colour = "black", alpha=1, size=2.5, angle = base_data$angle, inherit.aes = FALSE)

pdf("pltos.pdf", 6,6)
p + annotate("text", x=0, y=-min(scales)*4, label= "The number of publications\non machine learning to study\nthe effects of SNP", size = 2.5, fontface = "bold")
dev.off()
