library(ggplot2)
library(dplyr)

elbowK <- function(vet,
                          ylab = 'metric',
                          maxK = NULL,
                          main = 'Elbow method to find optimal k') {
#given any vector vet
#assuming it is always decreasing
#find analitically the elbow of the curve
  
  if (is.null(maxK)) {
    maxK <- length(vet)
  }
  Kopt <- 1:maxK
  data <- data.frame(K = Kopt, Metric = vet[Kopt])
  
  a <- data[1, ]
  b <- data[nrow(data), ]
  
  m <- -(abs(b$Metric - a$Metric)) / abs((b$K - a$K))
  q <- a$Metric - m * a$K
  mp <- -1 / m
  
  distances <- numeric(nrow(data))
  projection_points <- data.frame(K_proj = numeric(nrow(data)), Metric_proj = numeric(nrow(data)))
  
  for (i in 1:nrow(data)) {
    c_point <- data[i, ]
    qp <- c_point$Metric - mp * c_point$K
    K_proj <- (qp - q) / (m - mp)
    Metric_proj <- m * K_proj + q
    projection_points[i, ] <- c(K_proj, Metric_proj)
    distances[i] <- sqrt((c_point$K - K_proj)^2 + (c_point$Metric - Metric_proj)^2)
  }
  
  optimal_point <- data[which.max(distances), ]
  optimal_projection <- projection_points[which.max(distances), ]
  bestK <- optimal_point$K
  
    p <- ggplot(data, aes(x = K, y = Metric)) +
      geom_line(color = 'blue') +
      geom_point(color = 'blue') +
      annotate("segment", x = a$K, y = a$Metric, xend = b$K, yend = b$Metric,
               linetype = 'dashed', color = 'red') +
      annotate("segment", x = optimal_point$K, y = optimal_point$Metric,
               xend = optimal_projection$K_proj, yend = optimal_projection$Metric_proj,
               linetype = 'dashed', color = 'red') +
      geom_point(data = optimal_point, aes(x = K, y = Metric),
                 shape = 8, color = 'red', size = 3) +
      ggtitle(main) + ylab(ylab)+
      theme_bw()
    
    if (!is.null(maxK) && maxK > 0) {
      p <- p + xlim(0, maxK)
    }
    
    text_x <- quantile(data$K, 0.75)
    text_y <- max(data$Metric)
    p <- p + annotate("text", x = text_x, y = text_y,
                      label = paste('good k =', as.character(bestK)),
                      hjust = 0, vjust = 1)
    
    
    return(list(k=bestK, p=p, d=distances))
}
