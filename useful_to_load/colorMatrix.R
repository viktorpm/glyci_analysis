
colorMatrix = matrix(0,5,3)
rownames(colorMatrix) = c('blue','orange','green','purple','brown')
colnames(colorMatrix) = c('dark','midle','light')
colorMatrix['blue',] = c('#1D4871','#9EC2E9','#D6E9FC') 
colorMatrix['orange',] = c('#EB8104','#F5B887','#FFE3CC') 
colorMatrix['green',] = c('#127D49','#76EDB4','#CDFAE4') 
colorMatrix['purple',] = c('#AA67B2','#F2C3F7','#FCE3FF') 
colorMatrix['brown',] = c('#B66A34','#DBB79E','#FAE4D5') 

save(colorMatrix, file = 'f:/_R_WD/colorMatrix.RData')

