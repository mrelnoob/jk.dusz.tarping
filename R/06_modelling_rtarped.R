### Model diagnostics:
# One plot of residuals (https://www.r-bloggers.com/2011/07/model-validation-interpreting-residual-plots/):
plot(fitted(Cand.mod[[1]]), residuals(Cand.mod[[1]]),
     xlab = "Fitted Values", ylab = "Residuals")
abline(h=0, lty=2)
lines(smooth.spline(fitted(Cand.mod[[1]]), residuals(Cand.mod[[1]])))

# If mixed model:
# Check for residual pattern within groups (levels of random factor) and difference between groups
xyplot(residuals(glmm1) ~ fitted(glmm1) | Count$plot, main = "glmm1 – full model by plot",
       panel=function(x, y){
         panel.xyplot(x, y)
         panel.loess(x, y, span = 0.75)
         panel.lmline(x, y, lty = 2)  # Least squares broken line
       })


### EN cas de models mixtes, relire la vignette DHARMa, car ça change!
