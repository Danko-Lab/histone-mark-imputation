library(keras)

load("/local/workdir/zw355/proj/prj15-histone/models/H3k27me3.S5.V2.train.rdata.trainxy.rdata")

 
build_model <- function() {
  
  model <- keras_model_sequential() %>%
    layer_dense(units = 64, activation = "relu",
                input_shape = dim(x_train)[2]) %>%
    layer_dense(units = 64, activation = "relu") %>%
    layer_dense(units = 1)
  
  model %>% compile(
    loss = "mse",
    optimizer = optimizer_rmsprop(),
    metrics = list("mean_absolute_error")
  )
  
  model
}

model <- build_model()
model %>% summary()


# Display training progress by printing a single dot for each completed epoch.
print_dot_callback <- callback_lambda(
  on_epoch_end = function(epoch, logs) {
    if (epoch %% 5 == 0) cat("\n")
    cat(".")
  }
)    

epochs <- 500

# Fit the model and store training stats
history <- model %>% fit(
  x_train,
  y_train,
  epochs = epochs,
  validation_split = 0.2,
  verbose = 0,
  callbacks = list(print_dot_callback)
)

library(ggplot2)

plot(history, metrics = "mean_absolute_error", smooth = FALSE) +
  coord_cartesian(ylim = c(0, 5))
  
  early_stop <- callback_early_stopping(monitor = "val_loss", patience = 20)
  
  model <- build_model()
  history <- model %>% fit(
    x_train,
    y_train,
    epochs = epochs,
    validation_split = 0.2,
    verbose = 0,
    callbacks = list(early_stop, print_dot_callback)
  )
  
  plot(history, metrics = "mean_absolute_error", smooth = FALSE) +
  coord_cartesian(xlim = c(0, 150), ylim = c(0, 5))
  
  
c(loss, mae) %<-% (model %>% evaluate(x_test, y_test, verbose = 0))

paste0("Mean absolute error on test set: $", sprintf("%.2f", mae * 1000))

test_predictions <- model %>% predict(as.matrix(x_test))
test_predictions[ , 1]


