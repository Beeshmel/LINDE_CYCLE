import numpy as np
from CoolProp.CoolProp import PropsSI
from matplotlib import pyplot as plt
from numpy import arange
import CoolProp.CoolProp as CP

CP.set_config_string(CP.ALTERNATIVE_REFPROP_LIBRARY_PATH , 'C:\Program Files (x86)\REFPROP\REFPRP64.DLL')
print(CP.get_global_param_string("REFPROP_version"))
CP.set_config_bool(CP.REFPROP_USE_GERG , True)


def CompressorAdiabat(p0 , T0 , pk , st): #lk , Tk , Qk
    n_kol_stup_komp = 1
    while (pk / p0) ** (1 / n_kol_stup_komp) > st:
        n_kol_stup_komp += 1

    kpd_ad = 0.85
    R = 287.1  # газовая постоянная воздуха Дж/кг К
    nk = (pk / p0) ** (1 / n_kol_stup_komp)
 #   print("Степень сжатия в " , n_kol_stup_komp , "х-ступенчатом компрессоре" , round(nk , 2))
    Cp = PropsSI('C' , 'T' , T0 , 'P' , p0 , 'REFPROP::Air')
    Cv = PropsSI('CVMASS' , 'T' , T0 , 'P' , p0 , 'REFPROP::Air')
    n = Cp / Cv  # показатель адиабаты = показатель политропы 1.25 в некоторых источниках
    lk = n / (n - 1) * R * T0 * ((nk) ** ((n - 1) / n) - 1) / kpd_ad*n_kol_stup_komp  # политропное сжатие - рабочий процесс
    Tk = T0 * (1 + 1 / kpd_ad * ((nk) ** ((n - 1) / n) - 1))
    Qk = Cp * (Tk - T0)*n_kol_stup_komp
    print("Температура в конце сжатия " , round(Tk , 2) , "К")
    print("Удельная работа компрессора " , round(lk / 1000 , 2) , "кДж/кг")
    print("Теплота сжатия " , round(Qk / 1000 , 2) , "кДж/кг")
    return lk , Tk , Qk


def PumpIsoterm(p0 , T0 , pk ): #Lp , Tp
    k_pump = 0.75
    density = PropsSI("D" , "T" , T0 , "P|liquid" , (p0 + pk) / 2 , "REFPROP::Air")
    Lp = (pk - p0) / density / k_pump
    h0 = PropsSI("H" , "T" , T0 , "P|liquid" , p0 , "REFPROP::Air")
    hp = h0 + Lp
    Tp = PropsSI("T" , "H" , hp , "P" , pk , "REFPROP::Air")
    print("Температура после насоса " , round(Tp , 2) , "К")
    print("Работа насоса " , round(Lp / 1000 , 2) , " кДж/кг")
    return Lp , Tp


def Turbina(Toc , Qkm , Gtur , pt1 , Tt1 , pt2 , Tkm, st ):    # Ld_sum*k_tur
    n_st_det = 1
    while (pt1 / pt2) ** (1 / n_st_det) > st:
        n_st_det += 1

    k_tur = 0.75
    nd = (pt1 / pt2) ** (1 / n_st_det)
    print("Степень расширения " , nd , " Количество ступеней расширения " , n_st_det)
    Ld_sum = 0
    for i in arange(n_st_det):
        pt2 = pt1 / nd
        Cp = PropsSI('CPMASS' , 'T' , Tt1 , 'P' , (pt1+pt2)/2 , 'REFPROP::Air')
        T_heat = Qkm/Cp/Gtur+Tt1
        if T_heat > Tkm:
            T_heat = Tkm - 5
            if T_heat > 550:  # ограничение температуры по термомаслу
                T_heat = 550
        else:
            T_heat = Toc
        Cp = PropsSI('CPMASS' , 'T' , T_heat, 'P' , (pt1+pt2)/2 , 'REFPROP::Air')

        #print("Температура подогретого воздуха " , T_heat , " К")
        Qkm = Qkm - Cp * Gtur*(T_heat-Tt1)
        #print("Давление стало: " , pt2 / 1000000 , "МПа")
        Td , Ld = detander(pt1 , T_heat , pt2)
        Ld_sum = Ld_sum + Ld
        Tt1 = Td
        pt1 = pt2
    return Ld_sum*k_tur


 
