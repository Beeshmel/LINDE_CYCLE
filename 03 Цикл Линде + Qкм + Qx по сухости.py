from numpy import arange
from CoolProp.CoolProp import PropsSI
import mf
import CoolProp.CoolProp as CP
from matplotlib import pyplot as plt
import numpy as np

CP.set_config_string(CP.ALTERNATIVE_REFPROP_LIBRARY_PATH , 'C:\Program Files (x86)\REFPROP\REFPRP64.DLL')
print(CP.get_global_param_string("REFPROP_version"))
CP.set_config_bool(CP.REFPROP_USE_GERG , True)

R = 287.1  # газовая постоянная воздуха Дж/кг*К
k_mex = 0.9
dT = 5  # Недорекуперация на горячем конце принимается
NT = 10

T0 = 300
p0 = 101325  # 101325  # атмосферное давление

TIME_CHARGE = 8 * 60 * 60
TIME_DISCHARGE = 8 * 60 * 60
p7 = p8 = P_PUMP = 8 * 10 ** 6
Gmass = 100  # kg|s

masLx = []
masKPE = []
masPk = []
max_x = 0
KPE_MAX = 0
LXMIN = 10**10

T1 = T2 = T0
p1 = p4 = p5 = p6 = p0
h1 = PropsSI('H' , 'T' , T1 , 'P' , p1 , 'REFPROP::Air')  # Энтальпия прямого потока на входе в ТОА
s1 = PropsSI('S' , 'T' , T1 , 'P' , p1 , 'REFPROP::Air')

for pk in arange(15 * 10 ** 6 , 30 * 10 ** 6 , 1 * 10 ** 6):
    print("Давление" , pk/10 ** 6,'МПа')
    p2 = p3 = pk
    h2 = PropsSI('H' , 'T' , T2 , 'P' , p2 , 'REFPROP::Air')
    Lk , Tk , Qk = mf.CompressorAdiabat(p1 , T1 , p2 , 5)

    for T3 in arange(99 , 115 , 1):  # температура прямого потока на выходе
        h3 = PropsSI('H' , 'T' , T3 , 'P' , p3 , 'REFPROP::Air')

        h4 = h3
        T4 = PropsSI('T' , 'H' , h4 , 'P' , p4 , 'REFPROP::Air')

        Tf = T4
        hf = PropsSI('H' , 'T' , Tf , 'P|liquid' , p4 , 'REFPROP::Air')  # энтальпия жидкости после дросселя
        sf = PropsSI('S' , 'T' , Tf , 'P|liquid' , p4 , 'REFPROP::Air')

        T5 = Tf
        h5 = PropsSI('H' , 'T' , T5 , 'P|gas' , p5 , 'REFPROP::Air')

        x = 1 - PropsSI("Q","P",p4,"H",h4,"REFPROP::Air")

        T6 = T1 - dT
        h6 = PropsSI('H' , 'T' , T6 , 'P' , p6 , 'REFPROP::Air')

        Lp , T7 = mf.PumpIsoterm(p0 , Tf , p7)  # Температура после насоса
        h7 = h_pump = PropsSI('H' , 'T' , T7 , 'P' , p7 , 'Air')  # Энтальпия после насоса

        h8 = (h2-h3 - (h6-h5)*(1-x))/x +h7
        T8 = PropsSI('T' , 'H' , h8 , 'P' , p8 , 'REFPROP::Air')

        Zasechka = 0
        Qi_hot = (h2 - h3) / NT
        Qi_cold1 = (h6 - h5) / NT
        Qi_cold2 = (h8 - h7) / NT
        for i in arange(NT + 1):
            h2i = h2 - Qi_hot * i
            h6i = h6 - Qi_cold1 * i
            h8i = h8 - Qi_cold2 * i
            T2i = PropsSI('T' , 'H' , h2i , 'P' , p2 , 'REFPROP::Air')
            T6i = PropsSI('T' , 'H' , h6i , 'P' , p6 , 'REFPROP::Air')
            T8i = PropsSI('T' , 'H' , h8i , 'P' , p8 , 'REFPROP::Air')
            T68i = (T6i * (1 - x) + T8i * x)
            if (T2i < T68i)or(T2i < T6i)or(T2i < T8i):
                print("Пересечение в ТОА")
                Zasechka = 1
                break
            if (T2i - T68i < dT - 1):
                print("Недорекуперация в ТОА не достигнута")
                Zasechka = 1
                break

        if (Zasechka != 1):
            Gtur = x * Gmass
            Ltur = mf.Turbina(T0 - dT , 0 , Gtur , p8 , T8 , p0 , Tk, 5)  #наличие теплоты сжатия компрессора
            Wkm = Lk * Gmass
            Wp = Lp * Gtur / k_mex
            Wtur = Ltur * Gtur * k_mex
            KPE = (Wtur - Wp) / Wkm


            if KPE > KPE_MAX:

                KPE_MAX = KPE
                Qkm = Qk*Gmass
                Ltur_km = mf.Turbina(T0 - dT , Qkm , Gtur , p8 , T8 , p0 , Tk, 5)
                Wtur_km = Ltur_km * Gtur * k_mex
                KPE_km = (Wtur_km - Wp) / Wkm

                Lx = Lk / x

                max_x = x
                Eid = (h1 - hf) / (T1 * (s1 - sf) - (h1 - hf))
                qx = x*(h1 - hf)
                Ex = qx / Lk  # Холодильный коэффициент
                NUt = Ex / Eid  # Коэффициент термодинамического совершенства

                masPk.append(pk / 10 ** 6)
                masLx.append(Lx / 10 ** 3)
                masKPE.append(KPE)

                masRes = ["Температура T1=" , round(T1 , 2) , "K," , " давление " , p1 / 10 ** 6 , "МПа" ,
                          "Температура T2=" , round(T2 , 2) , "K," , " давление " , p2 / 10 ** 6 , "МПа" ,
                          "Температура T3=" , round(T3 , 2) , "K," , " давление " , p3 / 10 ** 6 , "МПа" ,
                          "Температура T4=" , round(T4 , 2) , "K," , " давление " , p4 / 10 ** 6 , "МПа" ,
                          "Температура T5=" , round(T5 , 2) , "K," , " давление " , p5 / 10 ** 6 , "МПа" ,
                          "Температура T6=" , round(T6 , 2) , "K," , " давление " , p6 / 10 ** 6 , "МПа",
                          "Температура T7=" , round(T7 , 2) , "K," , " давление " , p7 / 10 ** 6 , "МПа" ,
                          "Температура T8=" , round(T8 , 2) , "K," , " давление " , p8 / 10 ** 6 , "МПа"
                          ]


print(masRes)
print("Коэффициент ожижения " , round(max_x * 100 , 2) , "%")
print("Затраченная удельная работа " , (round(Lx / 10**3 , 3)) , " кДж/кг ж.в.")
print("Холодильный коэффициент " , round(Ex , 3))
print("Холодильный коэффициент ид. ож. цикла " , round(Eid , 3))
print("Коэффициент термодинамического совершенства " , round(NUt , 3))

print("1) Коэффициент преобразования энергии УАЭжв c учетом рекуперации " , KPE_MAX)
print("2) Коэффициент преобразования энергии УАЭжв c учетом рекуперации и теплоты сжатия КМ " , KPE_km)

x = np.array(masPk)
y = np.array(masLx)
plt.title("")
plt.xlabel("Давление нагнетания, МПа")
plt.ylabel("Удельная работа компрессора, кДж/кг ж.в.")
plt.plot(x , y , color="green")
plt.grid(True)
plt.show()

x = np.array(masPk)
y = np.array(masKPE)
plt.title("")
plt.xlabel("Давление нагнетания, МПа")
plt.ylabel("КПЭ")
plt.plot(x , y , color="green")
plt.grid(True)
plt.show()

# ['Температура T1=', 300, 'K,', ' давление ', 0.101325, 'МПа', 'Температура T2=', 300, 'K,', ' давление ', 22.0, 'МПа', 'Температура T3=', 103, 'K,', ' давление ', 22.0, 'МПа', 'Температура T4=', 79.66, 'K,', ' давление ', 0.101325, 'МПа', 'Температура T5=', 79.66, 'K,', ' давление ', 0.101325, 'МПа', 'Температура T6=', 295, 'K,', ' давление ', 0.101325, 'МПа', 'Температура T7=', 83.13, 'K,', ' давление ', 8.0, 'МПа', 'Температура T8=', 281.61, 'K,', ' давление ', 8.0, 'МПа']
# Коэффициент ожижения  71.27 %
# Затраченная удельная работа  933.304  кДж/кг ж.в.
# Холодильный коэффициент  0.455
# Холодильный коэффициент ид. ож. цикла  0.578
# Коэффициент термодинамического совершенства  0.788
# 1) Коэффициент преобразования энергии УАЭжв c учетом рекуперации  0.14676788348901054
# 2) Коэффициент преобразования энергии УАЭжв c учетом рекуперации и теплоты сжатия КМ  0.24349942940236255