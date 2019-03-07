from pysb import *
import pandas as pd
from pysb.core import *
from pysb.bng import *
from pysb.integrate import *
import matplotlib.pyplot as plt
import numpy as np
from scipy import interpolate
from scipy.interpolate import *
import scipy.interpolate
import numpy as np

# x = [0,1,23]
# y = [3798.81606,11838.08526, 5907.897248]
# plt.plot(x,y, marker ='o')
# plt.legend(loc = 'best')
# plt.show()

# df = pd.DataFrame({
# '1': [11838,
# 136357 ,
# 26292  ,
# 9109   ,
# 9346   ,
# 70389  ,
# 9228   ,
# 4328   ,
# 7432   ,
# 17032  ,
# 25557  ,
# 10675  ,
# 702    ,
# 4450   ,
# 2754   ,
# 5869],
# '23': [5907,
# 193589 ,
# 8123   ,
# 119922 ,
# 2808   ,
# 10537  ,
# 9943   ,
# 8173   ,
# 6773   ,
# 2290   ,
# 9845   ,
# 30566  ,
# 19469  ,
# 2812   ,
# 11913  ,
# 11197]}, index=['Casp8', 'Bid', 'Mcl1', 'Bax','Fadd',
# 'Bad',
# 'BclxL',
# 'Mlkl',
# 'Rip3',
# 'NFkB2',
# 'Tradd',
# 'Casp3',
# 'Casp7',
# 'IkBe',
# 'Nemo',
# 'RelA'])
#
# lines = df.plot.line(marker = 'o')
# plt.show()
# df1 = pd.read_csv("data_1_ms.csv")
# df2 = pd.read_csv("data_23_ms.csv")
df3 = pd.read_csv("ms_data_0124.csv")
# print(df1.head())
# print(df3['1'][0], df3['23'][0])
# quit()

# df1.plot.bar(x='1_prot', rot=0, title='1 hour', figsize=(22,10), color = ['red', 'blue'], fontsize=16)
# plt.savefig('1hr_cttr')
# df2.plot.bar(x='23_prot', rot=0, title='23 hour', figsize=(22,10),color = ['red', 'blue'], fontsize=16)
# # plt.savefig('23hr_cttr')

# df3.plot.bar(x='prot', rot=0, title='1 & 23 hr treated', figsize=(22,10),color = ['red', 'blue'],fontsize=16)
# plt.title('Label-Free MS: L929 cells treated with 10 ng/ml TNFa', fontsize = 18)
# plt.xlabel('Significant Proteins', fontsize=20)
# plt.ylabel('Amount [Molecules/Cell]', fontsize=20)
# plt.legend(loc = 'best',prop={'size': 20})
# plt.xticks(fontsize = 18)
# plt.yticks(fontsize = 18)
# plt.savefig('1_23hr_tr_hist')
# plt.show()
#
# # plt.plot([1,23],[11838, 5907], 'ro')
# # plt.xticks([1,23])
# # # plt.xlim(1, 23)
# # plt.legend('Casp8')
# # plt.show()
#
plt.figure(figsize=(20,10))
df3.plot.bar(x='1_prot', rot=0, title='Label-Free MS: L929 cells treated with 10 ng/ml TNFa', figsize=(22,10),color = ['red', 'green', 'blue'],fontsize=16)
# plt.scatter(df3['1_prot'][:], df3['CT'][:], label = '0 hr', color = 'red', lw =3)
# plt.scatter(df3['1_prot'][:], df3['1TR'][:], label = '1 hr', color = 'green', lw =3)
# plt.scatter(df3['1_prot'][:], df3['24TR'][:], label ='23 hr', color = 'blue',lw = 3)
# plt.scatter(x, pre_cal_kin, label ='preCal', color = 'green', lw = 1)
plt.xlabel('Significant Proteins', fontsize=20)
plt.ylabel('Amount [Molecules/Cell]', fontsize=20)
plt.xticks(fontsize = 18)
plt.yticks(fontsize = 18)
plt.title('Label-Free MS: L929 cells treated with 10 ng/ml TNFa', fontsize = 18)
plt.legend(loc = 'best',prop={'size': 20})
plt.savefig('ct_1_24hr_tr_pt_hist')
# plt.legend(flipnum, title = 'flip', loc=0, fontsize = 5)
plt.show()
#
# # plt.figure(figsize=(20,10))
# # # plt.scatter(x, params, label = 'PFu', color = 'blue', lw =1)
# # plt.scatter(df3['prot'][:], df3['norm1'][:], label = '1 hr', color = 'red',lw =3)
# # plt.scatter(df3['prot'][:], df3['norm23'][:], label ='23 hr', color = 'blue',lw = 3)
# # # plt.scatter(x, pre_cal_kin, label ='preCal', color = 'green', lw = 1)
# # plt.xlabel('Significant Proteins', fontsize=20)
# # plt.ylabel('Amount [Molecules/Cell]', fontsize=20)
# # plt.xticks(df3['prot'][:], fontsize = 18)
# # plt.yticks(fontsize = 18)
# # plt.title('Label-Free MS Normalized: L929 cells treated with 10 ng/ml TNFa', fontsize = 18)
# # plt.legend(loc = 'best',prop={'size': 20})
# # plt.savefig('1_23hr_tr_ptnorm')
# # # plt.legend(flipnum, title = 'flip', loc=0, fontsize = 5)
# # plt.show()
#
# # df3.set_index('prot').plot(subplots =True)
# # plt.show()
# #
# #
# # # X = df3[:,0]
# # # col_1= df3[:,1]
# # plt.plot(df3['prot'],df3['1'][:])
# # plt.show()

# fig = plt.figure()
# fig.subplots_adjust(hspace=0.4, wspace=0.4)
# # for i in range(1,7):
# #     plt.subplot(2,3,i)
# row = df3.iloc[0]
# row.plot(kind='scatter')
# # plt.scatter(df3.iloc[0])
# # plt.scatter(df3['23'][])
# # df3.set_index(['1','23'])
#     # plt.xticks(labels = ['1', '23'])
# plt.show()