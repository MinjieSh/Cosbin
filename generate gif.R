plot_data <- data
label <- lbl2
plot_data <-
  t(apply(plot_data, 1, function(x) x / sum(x)))

# ggtern(data = as.data.frame(plot_data),
#        aes(x=plot_data[,1], y=plot_data[,2], z=plot_data[,3], color = label)) +
#   geom_point(alpha = 0.3) +
#   labs( x       = "G1",
#         y       = "G2",
#         z       = "G3") + 
#   scale_color_manual(values=c("lightseagreen", "royalblue", "steelblue1", "green4", "red", "orange", "magenta")) +
#   theme(
#     title = element_text(size = 12, face = "bold"),
#     axis.title = element_text(size = 12, face = "bold"),
#     legend.text = element_text(size = 12),
#     strip.background = element_blank(),
#     strip.text.x = element_blank(),
#     strip.text.y = element_blank()
#   )


ggtern(data = as.data.frame(plot_data),
       aes(x=plot_data[,1], y=plot_data[,2], z=plot_data[,3], color = label)) +
  geom_point() +
  labs( x       = "G1",
        y       = "G2",
        z       = "G3") + 
  theme(
    title = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 12),
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    strip.text.y = element_blank()
  )


# Let's use the iris dataset
# iris
test <- data
#test <- data3
#test <- plot_data
#rs <- rowSums(test)
#test <- test[-which(rs > 1000), ]
#coltemp <- lbl2[-which(rs > 1000)]
coltemp <- lbl2

#coltemp[coltemp == 'sDEG G1'] <- "magenta"
#coltemp[coltemp == 'sDEG G2'] <- "orange"
#coltemp[coltemp == 'sDEG G3'] <- "red"
coltemp[coltemp == 'iCEG (labeled)'] <- "green4"
coltemp[coltemp == 'DEG (exclude sDEG)'] <- "steelblue1"
coltemp[coltemp == 'sDEG'] <- "magenta"
# Static chart
G1 <- test[,1]
G2 <- test[,2]
G3 <- test[,3]
#plot3d( G1, G2, G3)
plot3d( G1, G2, G3, col = coltemp, type = "p", r = 0.2)

plot3d( G1, G2, G3, col = coltemp, type = "p", r = 0.2,
        xlim=c(0,2000),
        ylim=c(0,2000),
        zlim=c(0,2000))

# We can indicate the axis and the rotation velocity
play3d( spin3d( axis = c(0, 0, 1), rpm = 20), duration = 10 )

# Save like gif
movie3d(
  movie="rotate", 
  spin3d( axis = c(0, 0, 1), rpm = 7),
  duration = 10, 
  type = "gif", 
  clean = TRUE,
  dir = file.path(getwd(), "3dplot")
)

#########################################

plot_data = data3
label = lbl2

i <- 5
plot_data <- outlis[[i]]

plot_data <-
  t(apply(plot_data, 1, function(x) x / sum(x)))

#emf(file ="tcctcc.emf", width = 6.4, height = 5)
ggtern(data = as.data.frame(plot_data),
       aes(x=plot_data[,1], y=plot_data[,2], z=plot_data[,3], color = label)) +
  geom_point(alpha = 0.3) +
  labs( x       = "G1",
        y       = "G2",
        z       = "G3") + 
  #scale_color_manual(values=c("green4", "steelblue1", "magenta", "orange", "red")) +
  scale_color_manual(values=c("magenta", "orange", "red", "green4", "steelblue1")) +
  theme(
    title = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 12),
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    strip.text.y = element_blank()
  )
# dev.off()

#################################################


###########################################################
dir_out <- file.path(getwd(), "images64")
dir.create(dir_out, recursive = TRUE)

result <- data3
row.names(result) <- 1:nrow(result)

#cos(arccos(1/sqrt(3)) - arccos(0.95))
threshold <- 0.728284

count <- 1
plot_lbl <- lbl2
while(max(cos_iDEG(result)) >= threshold){
  
  temp_SMG_ind <- which.max(cos_iDEG(result))
  result <- result[-temp_SMG_ind, ]
  result <- apply(result, 2, function(x) x / sum(x))
  
  plot_lbl <- plot_lbl[-temp_SMG_ind]
  
  if(count %% 100 == 1){
    plot_data <- result

    plot_data <-
      t(apply(plot_data, 1, function(x) x / sum(x)))

    p <-  ggtern(data = as.data.frame(plot_data),
                 aes(x=plot_data[,1], y=plot_data[,2], z=plot_data[,3], color = plot_lbl)) +
      geom_point(alpha = 0.3) +
      labs( x       = "G1",
            y       = "G2",
            z       = "G3") +
      scale_color_manual(values=c( "steelblue1","green4", "red", "red", "red")) +
      # scale_color_manual(values=c("magenta", "orange", "red", "green4", "steelblue1")) +
      theme(
        title = element_text(size = 12, face = "bold"),
        axis.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 12),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        strip.text.y = element_blank()
      )

    fp <- file.path(dir_out, paste0("0", (count-1) / 10, ".png"))
    ggsave(plot = p,
           filename = fp,
           device = "png",
           width = 8,
           height = 5
    )
    # print(p)
  }
  
  count <- count + 1
}

# print(dim(result))
# View(result)

# STEP 2 identify CEG, normalize based on CEG
# converge: # not change
count2 <- count + 1
while(length(which(cos_iCEG(result) >= 0.99)) != dim(result)[1]){
  temp_CEG_ind <- which(cos_iCEG(result) >= 0.99)
  result <- result[temp_CEG_ind,]
  result <- apply(result, 2, function(x) x / sum(x))
  # print(dim(result))
  plot_lbl <- plot_lbl[temp_CEG_ind]
  
  plot_data <- result
  
  plot_data <-
    t(apply(plot_data, 1, function(x) x / sum(x)))
  
  p <-  ggtern(data = as.data.frame(plot_data),
               aes(x=plot_data[,1], y=plot_data[,2], z=plot_data[,3], color = plot_lbl)) +
    geom_point(alpha = 0.3) +
    labs( x       = "G1",
          y       = "G2",
          z       = "G3") + 
    scale_color_manual(values=c("lightseagreen", "royalblue", "steelblue1", "green4", "red", "orange", "magenta")) +
    # scale_color_manual(values=c("magenta", "orange", "red", "green4", "steelblue1")) +
    theme(
      title = element_text(size = 12, face = "bold"),
      axis.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 12),
      strip.background = element_blank(),
      strip.text.x = element_blank(),
      strip.text.y = element_blank()
    )
  
  fp <- file.path(dir_out, paste0(count2, ".png"))
  ggsave(plot = p,
         filename = fp,
         device = "png",
         width = 8,
         height = 5
  )
  # print(p)
  
  count2 <- count2 + 1
  
  
}
print(dim(result))
# View(result)
data4 <- data3

# Final step normalization
ind <- as.numeric(row.names(result))
scalar <- colMeans(data4[ind, ] / result)
print(scalar)
for (i in 1:ncol(data4)) {
  data4[, i] <- data4[, i] / scalar[i]
}


plot_data <- data4

plot_data <-
  t(apply(plot_data, 1, function(x) x / sum(x)))

p <-  ggtern(data = as.data.frame(plot_data),
             aes(x=plot_data[,1], y=plot_data[,2], z=plot_data[,3], color = lbl2)) +
  geom_point(alpha = 0.3) +
  labs( x       = "G1",
        y       = "G2",
        z       = "G3") + 
  scale_color_manual(values=c("lightseagreen", "royalblue", "steelblue1", "green4", "red", "orange", "magenta")) +
  # scale_color_manual(values=c("magenta", "orange", "red", "green4", "steelblue1")) +
  theme(
    title = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 12),
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    strip.text.y = element_blank()
  )

fp <- file.path(dir_out, paste0(count2 + 1, ".png"))
ggsave(plot = p,
       filename = fp,
       device = "png",
       width = 8,
       height = 5
)
# print(p)


#####################################
library(magick)
# 
# imgs <- list.files(dir_out, full.names = TRUE)
# img_list <- lapply(imgs, image_read)
# 
# ## join the images together
# img_joined <- image_join(img_list)
# 
# ## animate at 2 frames per second
# img_animated <- image_animate(img_joined, fps = 4)
# 
# ## view animated image
# img_animated
# 
# ## save to disk
# image_write(image = img_animated,
#             path = "animation521asydirich.gif")

system("magick.exe convert -delay 40 *.png try.gif")