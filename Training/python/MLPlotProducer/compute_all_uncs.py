import statsmodels.stats.proportion as ssp

# Efficiencies
print("\n")
num = 17129.0
den = 63075
c_low, c_up = ssp.proportion_confint(num, den, alpha=1-0.68, method='beta')
print(("absolute efficiency after big or (both for evt and tau tuples)"))
print(("efficiency = {} \nvar_up = {}\nvar_down = {} \nunc_up = {}\nunc_down = {}").format(num/den, c_up, c_low, c_up-num/den,  num/den-c_low))
print("\n")
num =  15014.0
den =  63075
c_low, c_up = ssp.proportion_confint(num, den, alpha=1-0.68, method='beta')
print(("absolute efficiency after cut based"))
print(("efficiency = {} \nvar_up = {}\nvar_down = {} \nunc_up = {}\nunc_down = {}").format(num/den, c_up, c_low, c_up-num/den,  num/den-c_low))
print("\n")
num =  15014.0
den =  17129
c_low, c_up = ssp.proportion_confint(num, den, alpha=1-0.68, method='beta')
print(("algorithmic efficiency after cut based"))
print(("efficiency = {} \nvar_up = {}\nvar_down = {} \nunc_up = {}\nunc_down = {}").format(num/den, c_up, c_low, c_up-num/den,  num/den-c_low))
print("\n")
num =  15581
den =  63075
c_low, c_up = ssp.proportion_confint(num, den, alpha=1-0.68, method='beta')
print(("efficiency after DNN based"))
print(("efficiency = {} \nvar_up = {}\nvar_down = {} \nunc_up = {}\nunc_down = {}").format(num/den, c_up, c_low, c_up-num/den,  num/den-c_low))
print("\n")
num =  15581
den =  17129
c_low, c_up = ssp.proportion_confint(num, den, alpha=1-0.68, method='beta')
print(("algo efficiency after DNN based"))
print(("efficiency = {} \nvar_up = {}\nvar_down = {} \nunc_up = {}\nunc_down = {}").format(num/den, c_up, c_low, c_up-num/den,  num/den-c_low))

print("\n")

# Rates
num = 720582
den = 4016162
l1Tau_rate = 75818.34
BigOr_rate=  75818.34*720582/4016162

c_low, c_up = ssp.proportion_confint(num, den, alpha=1-0.68, method='beta')
print(("rate after big or (both for evt and tau tuples)"))
print(("efficiency = {} \nvar_up = {}\nvar_down = {} \nunc_up = {}\nunc_down = {}").format(num/den, c_up, c_low, c_up-num/den,  num/den-c_low))
print(("rate {} \nvar_up = {}\nvar_down = {} \nunc_up = {}\nunc_down = {}").format(num*l1Tau_rate/den,  c_up*l1Tau_rate, c_low*l1Tau_rate,  (c_up-num/den)*l1Tau_rate,(num/den-c_low)*l1Tau_rate))
print("\n")
num = 273923
den = 720582
c_low, c_up = ssp.proportion_confint(num, den, alpha=1-0.68, method='beta')
print(("rate after cut_based"))
print(("efficiency = {} \nvar_up = {}\nvar_down = {} \nunc_up = {}\nunc_down = {}").format(num/den, c_up, c_low, c_up-num/den,  num/den-c_low))
print(("rate = {} \nvar_up = {}\nvar_down = {} \nunc_up = {}\nunc_down = {}").format(num*BigOr_rate/den,  c_up*BigOr_rate, c_low*BigOr_rate,  (c_up-num/den)*BigOr_rate,(num/den-c_low)*BigOr_rate))
print("\n")

num = 273923
den = 4016162
c_low, c_up = ssp.proportion_confint(num, den, alpha=1-0.68, method='beta')
print(("absolute rate after cut_based"))
print(("efficiency = {} \nvar_up = {}\nvar_down = {} \nunc_up = {}\nunc_down = {}").format(num/den, c_up, c_low, c_up-num/den,  num/den-c_low))
print(("rate = {} \nvar_up = {}\nvar_down = {} \nunc_up = {}\nunc_down = {}").format(num*l1Tau_rate/den,  c_up*l1Tau_rate, c_low*l1Tau_rate,  (c_up-num/den)*l1Tau_rate,(num/den-c_low)*l1Tau_rate))
print("\n")

num = 264855
den = 720582
c_low, c_up = ssp.proportion_confint(num, den, alpha=1-0.68, method='beta')
print(("rate after dnn_score"))
print(("efficiency = {} \nvar_up = {}\nvar_down = {} \nunc_up = {}\nunc_down = {}").format(num/den, c_up, c_low, c_up-num/den,  num/den-c_low))
print(("rate {} \nvar_up = {}\nvar_down = {} \nunc_up = {}\nunc_down = {}").format(num*BigOr_rate/den,  c_up*BigOr_rate, c_low*BigOr_rate,  (c_up-num/den)*BigOr_rate,(num/den-c_low)*BigOr_rate))
print("\n")
num = 264855
den = 4016162
c_low, c_up = ssp.proportion_confint(num, den, alpha=1-0.68, method='beta')
print(("absolute rate after dnn_score"))
print(("efficiency = {} \nvar_up = {}\nvar_down = {} \nunc_up = {}\nunc_down = {}").format(num/den, c_up, c_low, c_up-num/den,  num/den-c_low))
print(("rate {} \nvar_up = {}\nvar_down = {} \nunc_up = {}\nunc_down = {}").format(num*l1Tau_rate/den,  c_up*l1Tau_rate, c_low*l1Tau_rate,  (c_up-num/den)*l1Tau_rate,(num/den-c_low)*l1Tau_rate))
