import pandas as pd
import rdkit
import signaturizer
import pandas as pd
import numpy as np
import scipy
from signaturizer import Signaturizer
from scipy.stats import pearsonr
import csv
from collections import defaultdict
with open(r'./finalDTP.csv', 'r') as file:
    reader = csv.reader(file)
    reader =pd.DataFrame(reader)


df=reader.drop(index=0,inplace=False)
groups = df.groupby([0])
# 按照 UniProt ID 分组，将每个 UniProt ID 对应的 smiles 转换为列表
drug_targets = dict(df.groupby([0])[2].apply(list))  ##0是smiles，2是uniprs
########################################################
import pandas as pd
import numpy as np
##import matplotlib.pyplot as plt
#from matplotlib.gridspec import GridSpec

ST_csv = pd.read_csv(r'./shared_all.csv')
nonST_csv = pd.read_csv(r'./non_all.csv')

ST_csv["TG"] = "ST"
nonST_csv["TG"] = "non-ST"

ST_allinone = pd.concat([ST_csv, nonST_csv])

### A1 to E5, P(s_i|ST), P(s_i|non-ST) and L(si) calculate
layers_we_consider = ['A1', 'A2', 'A3', 'A4', 'A5', 'B1', 'B2', 'B3', 'B4', 'B5',
                      'C1', 'C2', 'C3', 'C4', 'C5', 'D1', 'D2', 'D3', 'D4', 'D5',
                      'E1', 'E2', 'E3', 'E4', 'E5']
P_si_ST_all = pd.DataFrame()
P_si_nonST_all = pd.DataFrame()
L_si_all = pd.DataFrame()
#########
for layer in layers_we_consider:
    lower_limit = -1
    upper_limit = 1
    twenty_even_interval = np.linspace(lower_limit, upper_limit, num=20)
####共享
    P_si_ST = []
    for interval in range(1, 20):
        mask_ST = (ST_csv[layer] >= twenty_even_interval[interval-1]) & (ST_csv[layer] <= twenty_even_interval[interval])
        P_si_ST.append(mask_ST.sum() / ST_csv[layer].notna().sum())
    P_si_ST_all[layer] = P_si_ST
###非共享

    P_si_nonST = []
    for interval in range(1, 20):
        mask_nonST = (nonST_csv[layer] >= twenty_even_interval[interval-1]) & (nonST_csv[layer] <= twenty_even_interval[interval])
        P_si_nonST.append(mask_nonST.sum() / nonST_csv[layer].notna().sum())
    P_si_nonST_all[layer] = P_si_nonST
###似然比
    L_si = P_si_ST_all[layer] / P_si_nonST_all[layer]
    L_si_all[layer] = L_si



##全部的(没有print函数)
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
#import matplotlib.pyplot as plt
# 构建拟合后的指数函数

def exponential_func(x, a):
    return np.exp(a * x)
    
def fitted_exponential(x):
    return exponential_func(x, *popt)


fit_params = {}
for layer in layers_we_consider:
    x = twenty_even_interval[:-1]  # 去除最后一个元素作为自变量
    y = L_si_all[layer].values.flatten()  # 将似然比转换为一维数组作为因变量
    valid_indices = ~np.isnan(x) & ~np.isnan(y)
    similarity = x[valid_indices]
    likelihood = y[valid_indices]##################
    popt, pcov = curve_fit(exponential_func, similarity, likelihood)
    fit_params[layer] = popt
    x_fit = np.linspace(min(similarity), max(similarity), 100)
    y_fit = fitted_exponential(x_fit)


#########
data = pd.read_csv(r"./test.csv")
requirmentsmiles=list(data['smiles'])
#
r_sign_A1=Signaturizer('A1')
r_signatures_A1={}
r_sign_A2 = Signaturizer('A2')# ###2
r_signatures_A2= {}#1
r_sign_A3=Signaturizer('A3')
r_signatures_A3={}
r_sign_A4 = Signaturizer('A4')# ###2
r_signatures_A4 = {}#1
r_sign_A5 = Signaturizer('A5')# ###2
r_signatures_A5 = {}#1
######
r_sign_B1=Signaturizer('B1')
r_signatures_B1={}
r_sign_B2 = Signaturizer('B2')# ###2
r_signatures_B2= {}#1
r_sign_B3=Signaturizer('B3')
r_signatures_B3={}
r_sign_B4 = Signaturizer('B4')# ###2
r_signatures_B4 = {}#1
r_sign_B5 = Signaturizer('B5')# ###2
r_signatures_B5 = {}#1
######
r_sign_C1=Signaturizer('C1')
r_signatures_C1={}
r_sign_C2 = Signaturizer('C2')# ###2
r_signatures_C2= {}#1
r_sign_C3=Signaturizer('C3')
r_signatures_C3={}
r_sign_C4 = Signaturizer('C4')# ###2
r_signatures_C4 = {}#1
r_sign_C5 = Signaturizer('C5')# ###2
r_signatures_C5 = {}#1
###
r_sign_D1=Signaturizer('D1')
r_signatures_D1={}
r_sign_D2 = Signaturizer('D2')# ###2
r_signatures_D2= {}#1
r_sign_D3=Signaturizer('D3')
r_signatures_D3={}
r_sign_D4 = Signaturizer('D4')# ###2
r_signatures_D4 = {}#1
r_sign_D5 = Signaturizer('D5')# ###2
r_signatures_D5 = {}#1
##
r_sign_E1=Signaturizer('E1')
r_signatures_E1={}
r_sign_E2 = Signaturizer('E2')# ###2
r_signatures_E2= {}#1
r_sign_E3=Signaturizer('E3')
r_signatures_E3={}
r_sign_E4 = Signaturizer('E4')# ###2
r_signatures_E4 = {}#1
r_sign_E5 = Signaturizer('E5')# ###2
r_signatures_E5 = {}#1
for r_smi in requirmentsmiles:
    ##
    rsign_A1 = r_sign_A1.predict(r_smi)
    r_signatures_A1[r_smi] = rsign_A1.signature
    rsign_A2 = r_sign_A2.predict(r_smi)
    r_signatures_A2[r_smi] = rsign_A2.signature
    rsign_A3 = r_sign_A3.predict(r_smi)
    r_signatures_A3[r_smi] = rsign_A3.signature
    rsign_A4 = r_sign_A4.predict(r_smi)
    r_signatures_A4[r_smi] = rsign_A4.signature
    rsign_A5 = r_sign_A5.predict(r_smi)
    r_signatures_A5[r_smi] = rsign_A5.signature
##
    rsign_B1 = r_sign_B1.predict(r_smi)
    r_signatures_B1[r_smi] = rsign_B1.signature
    rsign_B2 = r_sign_B2.predict(r_smi)
    r_signatures_B2[r_smi] = rsign_B2.signature
    rsign_B3 = r_sign_B3.predict(r_smi)
    r_signatures_B3[r_smi] = rsign_B3.signature
    rsign_B4 = r_sign_B4.predict(r_smi)
    r_signatures_B4[r_smi] = rsign_B4.signature
    rsign_B5 = r_sign_B5.predict(r_smi)
    r_signatures_B5[r_smi] = rsign_B5.signature
    ###
    rsign_C1 = r_sign_C1.predict(r_smi)
    r_signatures_C1[r_smi] = rsign_C1.signature
    rsign_C2 = r_sign_C2.predict(r_smi)
    r_signatures_C2[r_smi] = rsign_C2.signature
    rsign_C3 = r_sign_C3.predict(r_smi)
    r_signatures_C3[r_smi] = rsign_C3.signature
    rsign_C4 = r_sign_C4.predict(r_smi)
    r_signatures_C4[r_smi] = rsign_C4.signature
    rsign_C5 = r_sign_C5.predict(r_smi)
    r_signatures_C5[r_smi] = rsign_C5.signature
    ##
    rsign_D1 = r_sign_D1.predict(r_smi)
    r_signatures_D1[r_smi] = rsign_D1.signature
    rsign_D2 = r_sign_D2.predict(r_smi)
    r_signatures_D2[r_smi] = rsign_D2.signature
    rsign_D3 = r_sign_D3.predict(r_smi)
    r_signatures_D3[r_smi] = rsign_D3.signature
    rsign_D4 = r_sign_D4.predict(r_smi)
    r_signatures_D4[r_smi] = rsign_D4.signature
    rsign_D5 = r_sign_D5.predict(r_smi)
    r_signatures_D5[r_smi] = rsign_D5.signature

    rsign_E1 = r_sign_E1.predict(r_smi)
    r_signatures_E1[r_smi] = rsign_E1.signature
    rsign_E2 = r_sign_E2.predict(r_smi)
    r_signatures_E2[r_smi] = rsign_E2.signature
    rsign_E3 = r_sign_E3.predict(r_smi)
    r_signatures_E3[r_smi] = rsign_E3.signature
    rsign_E4 = r_sign_E4.predict(r_smi)
    r_signatures_E4[r_smi] = rsign_E4.signature
    rsign_E5 = r_sign_E5.predict(r_smi)
    r_signatures_E5[r_smi] = rsign_E5.signature


sign_A1 = Signaturizer('A1')# ###2
sign_A2 = Signaturizer('A2')# ###2
sign_A3 = Signaturizer('A3')# ###2
sign_A4 = Signaturizer('A4')#1:sign_B4
sign_A5 = Signaturizer('A5')# ###2
sign_B1 = Signaturizer('B1')# ###2
sign_B2 = Signaturizer('B2')# ###2
sign_B3 = Signaturizer('B3')# ###2
sign_B4 = Signaturizer('B4')#1:sign_B4
sign_B5 = Signaturizer('B5')# ###2
sign_C1 = Signaturizer('C1')# ###2
sign_C2 = Signaturizer('C2')# ###2
sign_C3 = Signaturizer('C3')# ###2
sign_C4 = Signaturizer('C4')# ###2
sign_C5 = Signaturizer('C5')#1:sign_B4
sign_D1 = Signaturizer('D1')# ###2
sign_D2 = Signaturizer('D2')# ###2
sign_D3 = Signaturizer('D3')# ###2
sign_D4 = Signaturizer('D4')# ###2
sign_D5 = Signaturizer('D5')#1:sign_B4
sign_E1 = Signaturizer('E1')#1:sign_B4
sign_E2 = Signaturizer('E2')# ###2
sign_E3 = Signaturizer('E3')# ###2
sign_E4 = Signaturizer('E4')#1:sign_B4
sign_E5 = Signaturizer('E5')#1:sign_B4
##############
all_signatures_A1 = {}#1
all_signatures_A2 = {}#1
all_signatures_A3 = {}#1
all_signatures_A4 = {}#1
all_signatures_A5 = {}#1
all_signatures_B1 = {}#1
all_signatures_B2 = {}#1
all_signatures_B3 = {}#1
all_signatures_B4 = {}
all_signatures_B5 = {}#1
all_signatures_C1 = {}#1
all_signatures_C2 = {}#1
all_signatures_C3 = {}#1
all_signatures_C4 = {}#1
all_signatures_C5 = {}
all_signatures_D1 = {}
all_signatures_D2 = {}#1
all_signatures_D3 = {}
all_signatures_D4 = {}
all_signatures_D5 = {}
all_signatures_E1 = {}
all_signatures_E2 = {}#1
all_signatures_E3 = {}
all_signatures_E4 = {}
all_signatures_E5 = {}
for drug in drug_targets.keys():
    signature_A1 = sign_A1.predict(drug)#2：signature_B4 + sign_B4
    signature_A2 = sign_A2.predict(drug)#2：signature_B4 + sign_B4
    signature_A3 = sign_A3.predict(drug)#2：signature_B4 + sign_B4
    signature_A4 = sign_A4.predict(drug)#2：signature_B4 + sign_B4
    signature_A5 = sign_A5.predict(drug)#2：signature_B4 + sign_B4
    signature_B1 = sign_B1.predict(drug)#2：signature_B4 + sign_B4
    signature_B2 = sign_B2.predict(drug)#2：signature_B4 + sign_B4
    signature_B3 = sign_B3.predict(drug)#2：signature_B4 + sign_B4
    signature_B4 = sign_B4.predict(drug)#2：signature_B4 + sign_B4
    signature_B5 = sign_B5.predict(drug)#2：signature_B4 + sign_B4
    signature_C1 = sign_C1.predict(drug)#2：signature_B4 + sign_B4
    signature_C2 = sign_C2.predict(drug)#2：signature_B4 + sign_B4
    signature_C3 = sign_C3.predict(drug)#2：signature_B4 + sign_B4
    signature_C4 = sign_C4.predict(drug)#2：signature_B4 + sign_B4
    signature_C5 = sign_C5.predict(drug)#2：signature_B4 + sign_B4
    signature_D1 = sign_D1.predict(drug)#2：signature_B4 + sign_B4
    signature_D2 = sign_D2.predict(drug)#2：signature_B4 + sign_B4
    signature_D3 = sign_D3.predict(drug)#2：signature_B4 + sign_B4
    signature_D4 = sign_D4.predict(drug)#2：signature_B4 + sign_B4
    signature_D5 = sign_D5.predict(drug)#2：signature_B4 + sign_B4
    signature_E1 = sign_E1.predict(drug)#2：signature_B4 + sign_B4
    signature_E2 = sign_E2.predict(drug)#2：signature_B4 + sign_B4
    signature_E3 = sign_E3.predict(drug)#2：signature_B4 + sign_B4
    signature_E4 = sign_E4.predict(drug)#2：signature_B4 + sign_B4
    signature_E5 = sign_E5.predict(drug)#2：signature_B4 + sign_B4
    
    all_signatures_A1[drug] = signature_A1.signature#2 ：all_signatures_B4 + signature_B4
    all_signatures_A2[drug] = signature_A2.signature#2 ：all_signatures_B4 + signature_B4
    all_signatures_A3[drug] = signature_A3.signature#2 ：all_signatures_B4 + signature_B4
    all_signatures_A4[drug] = signature_A4.signature#2 ：all_signatures_B4 + signature_B4
    all_signatures_A5[drug] = signature_A5.signature#2 ：all_signatures_B4 + signature_B4
    all_signatures_B1[drug] = signature_B1.signature#2 ：all_signatures_B4 + signature_B4
    all_signatures_B2[drug] = signature_B2.signature#2 ：all_signatures_B4 + signature_B4
    all_signatures_B3[drug] = signature_B3.signature#2 ：all_signatures_B4 + signature_B4
    all_signatures_B4[drug] = signature_B4.signature#2 ：all_signatures_B4 + signature_B4
    all_signatures_B5[drug] = signature_B5.signature#2 ：all_signatures_B4 + signature_B4
    all_signatures_C1[drug] = signature_C1.signature#2 ：all_signatures_B4 + signature_B4
    all_signatures_C2[drug] = signature_C2.signature#2 ：all_signatures_B4 + signature_B4
    all_signatures_C3[drug] = signature_C3.signature#2 ：all_signatures_B4 + signature_B4
    all_signatures_C4[drug] = signature_C4.signature#2 ：all_signatures_B4 + signature_B4
    all_signatures_C5[drug] = signature_C5.signature#2 ：all_signatures_B4 + signature_B4
    all_signatures_D1[drug] = signature_D1.signature#2 ：all_signatures_B4 + signature_B4
    all_signatures_D2[drug] = signature_D2.signature#2 ：all_signatures_B4 + signature_B4
    all_signatures_D3[drug] = signature_D3.signature#2 ：all_signatures_B4 + signature_B4
    all_signatures_D4[drug] = signature_D4.signature#2 ：all_signatures_B4 + signature_B4
    all_signatures_D5[drug] = signature_D5.signature#2 ：all_signatures_B4 + signature_B4
    all_signatures_E1[drug] = signature_E1.signature#2 ：all_signatures_B4 + signature_B4
    all_signatures_E2[drug] = signature_E2.signature#2 ：all_signatures_B4 + signature_B4
    all_signatures_E3[drug] = signature_E3.signature#2 ：all_signatures_B4 + signature_B4
    all_signatures_E4[drug] = signature_E4.signature#2 ：all_signatures_B4 + signature_B4
    all_signatures_E5[drug] = signature_E5.signature#2 ：all_signatures_B4 + signature_B4

# import matplotlib.pyplot as plt
# 创建一个字典来存储不同 cutoff 下的 hit_rate
cutoff = 1000
best_T = {}
for i in drug_targets3.keys():
    signature1_A1 = r_signatures_A1[i]#2
    signature1_A2 = r_signatures_A2[i]#2
    signature1_A3 = r_signatures_A3[i]#2
    signature1_A4 = r_signatures_A4[i]#2
    signature1_A5 = r_signatures_A5[i]#2

    signature1_B1 = r_signatures_B1[i]#2
    signature1_B2 = r_signatures_B2[i]#2
    signature1_B3 = r_signatures_B3[i]#2
    signature1_B4 = r_signatures_B4[i]#2
    signature1_B5 = r_signatures_B5[i]#2

    signature1_C1 = r_signatures_C1[i]#2
    signature1_C2 = r_signatures_C2[i]#2
    signature1_C3 = r_signatures_C3[i]#2
    signature1_C4 = r_signatures_C4[i]#2
    signature1_C5 = r_signatures_C5[i]#2

    signature1_D1 = r_signatures_D1[i]#2
    signature1_D2 = r_signatures_D2[i]#2
    signature1_D3 = r_signatures_D3[i]#2
    signature1_D4 = r_signatures_D4[i]#2
    signature1_D5 = r_signatures_D5[i]#2

    signature1_E1 = r_signatures_E1[i]#2
    signature1_E2 = r_signatures_E2[i]#2
    signature1_E3 = r_signatures_E3[i]#2
    signature1_E4 = r_signatures_E4[i]#2
    signature1_E5 = r_signatures_E5[i]#2

    lsis=defaultdict(list)####
    for j,jtarget in drug_targets.items():
        signature2_A1 = all_signatures_A1[j]#2
        signature2_A2 = all_signatures_A2[j]#2
        signature2_A3 = all_signatures_A3[j]#2
        signature2_A4 = all_signatures_A4[j]#2
        signature2_A5 = all_signatures_A5[j]#2

        signature2_B1 = all_signatures_B1[j]#2
        signature2_B2 = all_signatures_B2[j]#2
        signature2_B3 = all_signatures_B3[j]#2
        signature2_B4 = all_signatures_B4[j]#2
        signature2_B5 = all_signatures_B5[j]#2

        signature2_C1 = all_signatures_C1[j]#2
        signature2_C2 = all_signatures_C2[j]#2
        signature2_C3 = all_signatures_C3[j]#2
        signature2_C4 = all_signatures_C4[j]#2
        signature2_C5 = all_signatures_C5[j]#2

        signature2_D1 = all_signatures_D1[j]#2
        signature2_D2 = all_signatures_D2[j]#2
        signature2_D3 = all_signatures_D3[j]#2
        signature2_D4 = all_signatures_D4[j]#2
        signature2_D5 = all_signatures_D5[j]#2

        signature2_E1 = all_signatures_E1[j]#2
        signature2_E2 = all_signatures_E2[j]#2
        signature2_E3 = all_signatures_E3[j]#2
        signature2_E4 = all_signatures_E4[j]#2
        signature2_E5 = all_signatures_E5[j]#2

        if j==i:
            continue
        correlation_A1 = np.corrcoef(signature1_A1, signature2_A1)[0, 1]##3
        correlation_A2 = np.corrcoef(signature1_A2, signature2_A2)[0, 1]
        correlation_A3 = np.corrcoef(signature1_A3, signature2_A3)[0, 1]##3
        correlation_A4 = np.corrcoef(signature1_A4, signature2_A4)[0, 1]
        correlation_A5 = np.corrcoef(signature1_A5, signature2_A5)[0, 1]##3

        correlation_B1 = np.corrcoef(signature1_B1, signature2_B1)[0, 1]
        correlation_B2 = np.corrcoef(signature1_B2, signature2_B2)[0, 1]##3
        correlation_B3 = np.corrcoef(signature1_B3, signature2_B3)[0, 1]
        correlation_B4 = np.corrcoef(signature1_B4, signature2_B4)[0, 1]##3
        correlation_B5 = np.corrcoef(signature1_B5, signature2_B5)[0, 1]

        correlation_C1 = np.corrcoef(signature1_C1, signature2_C1)[0, 1]
        correlation_C2 = np.corrcoef(signature1_C2, signature2_C2)[0, 1]##3
        correlation_C3 = np.corrcoef(signature1_C3, signature2_C3)[0, 1]
        correlation_C4 = np.corrcoef(signature1_C4, signature2_C4)[0, 1]##3
        correlation_C5 = np.corrcoef(signature1_C5, signature2_C5)[0, 1]

        correlation_D1 = np.corrcoef(signature1_D1, signature2_D1)[0, 1]
        correlation_D2 = np.corrcoef(signature1_D2, signature2_D2)[0, 1]##3
        correlation_D3 = np.corrcoef(signature1_D3, signature2_D3)[0, 1]
        correlation_D4 = np.corrcoef(signature1_D4, signature2_D4)[0, 1]##3
        correlation_D5 = np.corrcoef(signature1_D5, signature2_D5)[0, 1]

        correlation_E1 = np.corrcoef(signature1_E1, signature2_E1)[0, 1]
        correlation_E2 = np.corrcoef(signature1_E2, signature2_E2)[0, 1]##3
        correlation_E3 = np.corrcoef(signature1_E3, signature2_E3)[0, 1]
        correlation_E4 = np.corrcoef(signature1_E4, signature2_E4)[0, 1]##3
        correlation_E5 = np.corrcoef(signature1_E5, signature2_E5)[0, 1]

        layer_to_predict = "A1"  # 要预测的层#########@@@@@@@@改8
        if layer_to_predict in fit_params:
            popt = fit_params[layer_to_predict]  # 获取拟合参数
            lsiA1= fitted_exponential(correlation_A1)  # 使用新的 x 值计算新的 y 值###########@@@@@@@@@改9，10
        layer_to_predict = "A2"  # 要预测的层#########@@@@@@@@改8
        if layer_to_predict in fit_params:
            popt = fit_params[layer_to_predict]  # 获取拟合参数
            lsiA2 = fitted_exponential(correlation_A2)  # 使用新的 x 值计算新的 y 值###########@@@@@@@@@改9，10
        layer_to_predict = "A3"  # 要预测的层#########@@@@@@@@改8
        if layer_to_predict in fit_params:
            popt = fit_params[layer_to_predict]  # 获取拟合参数
            lsiA3 = fitted_exponential(correlation_A3)  # 使用新的 x 值计算新的 y 值###########@@@@@@@@@改9，10
        layer_to_predict = "A4"  # 要预测的层#########@@@@@@@@改8
        if layer_to_predict in fit_params:
            popt = fit_params[layer_to_predict]  # 获取拟合参数
            lsiA4 = fitted_exponential(correlation_A4)  # 使用新的 x 值计算新的 y 值###########@@@@@@@@@改9，10
        layer_to_predict = "A5"  # 要预测的层#########@@@@@@@@改8
        if layer_to_predict in fit_params:
            popt = fit_params[layer_to_predict]  # 获取拟合参数
            lsiA5 = fitted_exponential(correlation_A5)  # 使用新的 x 值计算新的 y 值###########@@@@@@@@@改9，10

        layer_to_predict = "B1"  # 要预测的层#########@@@@@@@@改8
        if layer_to_predict in fit_params:
            popt = fit_params[layer_to_predict]  # 获取拟合参数
            lsiB1= fitted_exponential(correlation_B1)  # 使用新的 x 值计算新的 y 值###########@@@@@@@@@改9，10
        layer_to_predict = "B2"  # 要预测的层#########@@@@@@@@改8
        if layer_to_predict in fit_params:
            popt = fit_params[layer_to_predict]  # 获取拟合参数
            lsiB2 = fitted_exponential(correlation_B2)  # 使用新的 x 值计算新的 y 值###########@@@@@@@@@改9，10
        layer_to_predict = "B3"  # 要预测的层#########@@@@@@@@改8
        if layer_to_predict in fit_params:
            popt = fit_params[layer_to_predict]  # 获取拟合参数
            lsiB3 = fitted_exponential(correlation_B3)  # 使用新的 x 值计算新的 y 值###########@@@@@@@@@改9，10
        layer_to_predict = "B4"  # 要预测的层#########@@@@@@@@改8
        if layer_to_predict in fit_params:
            popt = fit_params[layer_to_predict]  # 获取拟合参数
            lsiB4 = fitted_exponential(correlation_B4)  # 使用新的 x 值计算新的 y 值###########@@@@@@@@@改9，10
        layer_to_predict = "B5"  # 要预测的层#########@@@@@@@@改8
        if layer_to_predict in fit_params:
            popt = fit_params[layer_to_predict]  # 获取拟合参数
            lsiB5 = fitted_exponential(correlation_B5)  # 使用新的 x 值计算新的 y 值###########@@@@@@@@@改9，10

        layer_to_predict = "C1"  # 要预测的层#########@@@@@@@@改8
        if layer_to_predict in fit_params:
            popt = fit_params[layer_to_predict]  # 获取拟合参数
            lsiC1= fitted_exponential(correlation_C1)  # 使用新的 x 值计算新的 y 值###########@@@@@@@@@改9，10
        layer_to_predict = "C2"  # 要预测的层#########@@@@@@@@改8
        if layer_to_predict in fit_params:
            popt = fit_params[layer_to_predict]  # 获取拟合参数
            lsiC2 = fitted_exponential(correlation_C2)  # 使用新的 x 值计算新的 y 值###########@@@@@@@@@改9，10
        layer_to_predict = "C3"  # 要预测的层#########@@@@@@@@改8
        if layer_to_predict in fit_params:
            popt = fit_params[layer_to_predict]  # 获取拟合参数
            lsiC3 = fitted_exponential(correlation_C3)  # 使用新的 x 值计算新的 y 值###########@@@@@@@@@改9，10
        layer_to_predict = "C4"  # 要预测的层#########@@@@@@@@改8
        if layer_to_predict in fit_params:
            popt = fit_params[layer_to_predict]  # 获取拟合参数
            lsiC4 = fitted_exponential(correlation_C4)  # 使用新的 x 值计算新的 y 值###########@@@@@@@@@改9，10
        layer_to_predict = "C5"  # 要预测的层#########@@@@@@@@改8
        if layer_to_predict in fit_params:
            popt = fit_params[layer_to_predict]  # 获取拟合参数
            lsiC5 = fitted_exponential(correlation_C5)  # 使用新的 x 值计算新的 y 值###########@@@@@@@@@改9，10

        layer_to_predict = "D1"  # 要预测的层#########@@@@@@@@改8
        if layer_to_predict in fit_params:
            popt = fit_params[layer_to_predict]  # 获取拟合参数
            lsiD1= fitted_exponential(correlation_D1)  # 使用新的 x 值计算新的 y 值###########@@@@@@@@@改9，10
        layer_to_predict = "D2"  # 要预测的层#########@@@@@@@@改8
        if layer_to_predict in fit_params:
            popt = fit_params[layer_to_predict]  # 获取拟合参数
            lsiD2 = fitted_exponential(correlation_D2)  # 使用新的 x 值计算新的 y 值###########@@@@@@@@@改9，10
        layer_to_predict = "D3"  # 要预测的层#########@@@@@@@@改8
        if layer_to_predict in fit_params:
            popt = fit_params[layer_to_predict]  # 获取拟合参数
            lsiD3 = fitted_exponential(correlation_D3)  # 使用新的 x 值计算新的 y 值###########@@@@@@@@@改9，10
        layer_to_predict = "D4"  # 要预测的层#########@@@@@@@@改8
        if layer_to_predict in fit_params:
            popt = fit_params[layer_to_predict]  # 获取拟合参数
            lsiD4 = fitted_exponential(correlation_D4)  # 使用新的 x 值计算新的 y 值###########@@@@@@@@@改9，10
        layer_to_predict = "D5"  # 要预测的层#########@@@@@@@@改8
        if layer_to_predict in fit_params:
            popt = fit_params[layer_to_predict]  # 获取拟合参数
            lsiD5 = fitted_exponential(correlation_D5)  # 使用新的 x 值计算新的 y 值###########@@@@@@@@@改9，10

        layer_to_predict = "E1"  # 要预测的层#########@@@@@@@@改8
        if layer_to_predict in fit_params:
            popt = fit_params[layer_to_predict]  # 获取拟合参数
            lsiE1= fitted_exponential(correlation_E1)  # 使用新的 x 值计算新的 y 值###########@@@@@@@@@改9，10
        layer_to_predict = "E2"  # 要预测的层#########@@@@@@@@改8
        if layer_to_predict in fit_params:
            popt = fit_params[layer_to_predict]  # 获取拟合参数
            lsiE2 = fitted_exponential(correlation_E2)  # 使用新的 x 值计算新的 y 值###########@@@@@@@@@改9，10
        layer_to_predict = "E3"  # 要预测的层#########@@@@@@@@改8
        if layer_to_predict in fit_params:
            popt = fit_params[layer_to_predict]  # 获取拟合参数
            lsiE3 = fitted_exponential(correlation_E3)  # 使用新的 x 值计算新的 y 值###########@@@@@@@@@改9，10
        layer_to_predict = "E4"  # 要预测的层#########@@@@@@@@改8
        if layer_to_predict in fit_params:
            popt = fit_params[layer_to_predict]  # 获取拟合参数
            lsiE4 = fitted_exponential(correlation_E4)  # 使用新的 x 值计算新的 y 值###########@@@@@@@@@改9，10
        layer_to_predict = "E5"  # 要预测的层#########@@@@@@@@改8
        if layer_to_predict in fit_params:
            popt = fit_params[layer_to_predict]  # 获取拟合参数
            lsiE5 = fitted_exponential(correlation_E5)  # 使用新的 x 值计算新的 y 值###########@@@@@@@@@改9，10
        lsi=lsiA1*lsiA2*lsiA3*lsiA4*lsiA5*lsiB1*lsiB2*lsiB3*lsiB4*lsiB5*lsiC1*lsiC2*lsiC3*lsiC4*lsiC5*lsiD1*lsiD2*lsiD3*lsiD4*lsiD5*lsiE1*lsiE2*lsiE3*lsiE4*lsiE5######
        if lsi > cutoff:
            for eachtarget in jtarget:
                lsis[eachtarget].append(lsi)

####@@
    target_scores = {}
    for T, lsi_scores in lsis.items():
        avg_score = np.mean(lsi_scores)  # 计算平均分数
        target_scores[T] = avg_score
        

# 按照平均分数进行排序
    sorted_targets = sorted(target_scores.items(), key=lambda x: x[1], reverse=True)

# 选择前三个靶点作为 best_target
    best_target_candidates = [target for target, score in sorted_targets[:1]]
    best_T[i] = best_target_candidates

best_target_candidates


