# https://medium.com/@andriylazorenko/tensorflow-performance-test-cpu-vs-gpu-79fcd39170c

# https://towardsdatascience.com/keras-with-r-predicting-car-sales-31f48a58bf6


#export PATH=/programs/miniconda3/bin:$PATH
#source activate tensorflow-gpu
#
#python
#
#>>import tensorflow as tf
#
#source deactivate


library(keras)

setwd("/home/michaeljgrogan/Documents/a_documents/computing/data science/datasets")

cars<-read.csv("cars.csv")

normalize <- function(x) {
  return ((x - min(x)) / (max(x) - min(x)))
}

maxmindf <- as.data.frame(lapply(cars, normalize))
attach(maxmindf)

train_index <- sample(1:nrow(maxmindf), 0.8 * nrow(maxmindf))
test_index <- setdiff(1:nrow(maxmindf), train_index)

# Build X_train, y_train, X_test, y_test
X_train <- as.matrix(maxmindf[train_index, -15])
y_train <- as.matrix(maxmindf[train_index, "sales"])

X_test <- as.matrix(maxmindf[test_index, -15])
y_test <- as.matrix(maxmindf[test_index, "sales"])

model <- keras_model_sequential() 

model %>% 
  layer_dense(units = 12, activation = 'relu', kernel_initializer='RandomNormal', input_shape = c(6)) %>% 
  layer_dense(units = 8, activation = 'relu') %>%
  layer_dense(units = 1, activation = 'linear')

summary(model)

model %>% compile(
  loss = 'mean_squared_error',
  optimizer = 'adam',
  metrics = c('mae')
)

history <- model %>% fit(
  X_train, y_train, 
  epochs = 150, batch_size = 50, 
  validation_split = 0.2
)

model %>% evaluate(X_test, y_test)

pred <- data.frame(y = predict(model, as.matrix(X_test)))
df<-data.frame(pred,X_test)
attach(df)
deviation=((pred-sales)/sales)
mean(deviation$y)*100
