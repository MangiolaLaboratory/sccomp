# setup inverse_multinomial_logistic
inverse_multinomial_logistic = c(-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5)

# transform to probability
probability_vector = softmax(inverse_multinomial_logistic )

# Plot
plot(
  inverse_multinomial_logistic ,
  boot::logit(probability_vector)
)

# Regression
lm( boot::logit(probability_vector) ~ inverse_multinomial_logistic )

# Add line
abline(-5.331, 1.058  )
