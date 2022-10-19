###Fouling vs fouled number
#making table
observed_table <- matrix(c(423, 82, 46, 11), nrow = 2, ncol = 2, byrow = T)
rownames(observed_table) <- c('QM.fouled', 'ZM.fouled')
colnames(observed_table) <- c('QM.fouling', 'ZM.fouling')
observed_table

#chi-sq test
chi <- chisq.test(observed_table)
chi


observed_table2 <- matrix(c(12, 43, 47, 317), nrow = 2, ncol = 2, byrow = T)
rownames(observed_table2) <- c('QM.fouled', 'ZM.fouled')
colnames(observed_table2) <- c('QM.fouling', 'ZM.fouling')
observed_table2
chi2 <- chisq.test(observed_table2)
chi2

observed_table3 <- matrix(c(435, 125, 93, 328), nrow = 2, ncol = 2, byrow = T)
rownames(observed_table3) <- c('QM.fouled', 'ZM.fouled')
colnames(observed_table3) <- c('QM.fouling', 'ZM.fouling')
observed_table3
chi3 <- chisq.test(observed_table3)
chi3
chi3$residuals

require(corrplot)
corrplot(chi3$residuals, is.cor = FALSE)
