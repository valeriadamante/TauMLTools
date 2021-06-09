import statsmodels.stats.proportion as ssp

l1TauRate = 75818.34
bigOrRate=  75818.34*720582/4016162
lumiFraction = 2./1.68


def EvaluateRate(num, den, doL1=False, doBOR=False, doLumiF=False):
    l1TauRate = 75818.34
    bigOrRate=  75818.34*720582/4016162
    lumiFraction = 2./1.68
    eff = num/den
    c_low, c_up = ssp.proportion_confint(num, den, alpha=1-0.68, method='beta')
    print(("total numerator = \t {} \ntotal denominator = \t {} \nefficiency = \t {} \nvar_up = \t {} \nvar_down = \t {} \nunc_up = \t {} \nunc_down = \t {}").format(num, den, eff, c_up, c_low, c_up-eff, eff-c_low))
    if(doL1==True):
        eff_rate = eff*l1TauRate
        c_low_rate = c_low*l1TauRate
        c_up_rate = c_up*l1TauRate
        print(("\nL1 Rate = \t {} \nvar_up = \t {} \nvar_down = \t {} \nunc_up \t {} \nunc_down = \t {}").format(eff_rate,  c_up_rate, c_low_rate, c_up_rate-eff_rate, eff_rate-c_low_rate))
    if(doBOR==True):
        eff_rate = eff*bigOrRate
        c_low_rate = c_low*bigOrRate
        c_up_rate = c_up*bigOrRate
        print(("\nBigOR Rate = \t {} \nvar_up = \t {} \nvar_down = \t {} \nunc_up \t {} \nunc_down = \t {}").format(eff_rate,  c_up_rate, c_low_rate, c_up_rate-eff_rate, eff_rate-c_low_rate))
    if(doLumiF==True):
        eff_lumi = eff*lumiFraction
        c_low_lumi = c_low*lumiFraction
        c_up_lumi = c_up*lumiFraction
        print(("\nLumiFraction = \t {} \nvar_up = \t {} \nvar_down = \t {} \nunc_up \t {} \nunc_down = \t {}").format(eff_lumi,  c_up_lumi, c_low_lumi, c_up_lumi-eff_lumi, eff_lumi-c_low_lumi))
    if(doLumiF==True and doL1==True):
        eff_rate = eff*lumiFraction*l1TauRate
        c_low_rate = c_low*lumiFraction*l1TauRate
        c_up_rate = c_up*lumiFraction*l1TauRate
        print(("\nLumiFraction & L1 Rate = \t {} \nvar_up = \t {} \nvar_down = \t {} \nunc_up \t {} \nunc_down = \t {}").format(eff_rate,  c_up_rate, c_low_rate, c_up_rate-eff_rate, eff_rate-c_low_rate))
    if(doLumiF==True and doBOR==True):
        eff_rate = eff*lumiFraction*bigOrRate
        c_low_rate = c_low*lumiFraction*bigOrRate
        c_up_rate = c_up*lumiFraction*bigOrRate
        print(("\nLumiFraction & BigOR Rate = \t {} \nvar_up = \t {} \nvar_down = \t {} \nunc_up \t {} \nunc_down = \t {}").format(eff_rate,  c_up_rate, c_low_rate, c_up_rate-eff_rate, eff_rate-c_low_rate))


# RATES

print(("bigOrRate = {}, L1TauRate = {}, LumiFraction = {}").format(bigOrRate, l1TauRate, lumiFraction))

# cutBased
print("cutBased Rate")
num = 273923
den = 720582

EvaluateRate(num, den, False, True, True)
print("\n")
den = 4016162
EvaluateRate(num, den, True, False, True)
print("\n")

# 5 kHz
print("5 kHz Rate ")
num = 264855
den = 720582
EvaluateRate(num, den, False, True, True)
print("\n")
den = 4016162
EvaluateRate(num, den, True, False, True)
print("\n")


# 4 kHz
print("4 kHz Rate ")
num = 211884
den = 720582
EvaluateRate(num, den, False, True, True)
print("\n")
den = 4016162
EvaluateRate(num, den, True, False, True)
print("\n")

# 3 kHz
print("3 kHz Rate ")
num = 158913
den = 720582
EvaluateRate(num, den, False, True, True)
print("\n")
den = 4016162
EvaluateRate(num, den, True, False, True)
print("\n")

'''



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


'''
