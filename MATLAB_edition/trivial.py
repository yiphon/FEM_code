from xlrd import *

a=open_workbook('axisy_eg_input.xlsx')
source=a.sheets()[0]
nrow=source.nrows
p=[]

for i in range(nrow):
    p.append(source.row_values(i))
    print(source.row_values(i))
print(p)