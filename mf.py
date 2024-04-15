import numpy as np
from CoolProp.CoolProp import PropsSI
from matplotlib import pyplot as plt
from numpy import arange
import CoolProp.CoolProp as CP

CP.set_config_string(CP.ALTERNATIVE_REFPROP_LIBRARY_PATH , 'C:\Program Files (x86)\REFPROP\REFPRP64.DLL')
print(CP.get_global_param_string("REFPROP_version"))
CP.set_config_bool(CP.REFPROP_USE_GERG , True)

def KolStupDet(p0 , pk , st):
    n_kol_stup_det = 1
    while (p0 / pk) ** (1 / n_kol_stup_det) > st:
        n_kol_stup_det += 1
    # print("Количество ступеней детандера" , n_kol_stup_det)
    return n_kol_stup_det

def KolStupKomp(p0 , pk , st): #лишняяя функция УДАЛИТЬ!!!
    n_kol_stup_komp = 1
    while (pk / p0) ** (1 / n_kol_stup_komp) > st:
        n_kol_stup_komp += 1
    # print("Количество ступеней компрессора" , n_kol_stup_komp)
    return n_kol_stup_komp

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

def Compressor(p0 , T0 , pk,composition,st): #lk , Tk , Qk
    n_kol_stup_komp = 1
    while (pk / p0) ** (1 / n_kol_stup_komp) > st:
        n_kol_stup_komp += 1
    kpd_ad = 0.85
    R = 287.1  # газовая постоянная воздуха Дж/кг К
    nk = (pk / p0) ** (1 / n_kol_stup_komp)
    print("Степень сжатия в " , n_kol_stup_komp , "х-ступенчатом компрессоре" , round(nk , 2))
    Cp = PropsSI('C' , 'T' , T0 , 'P' , p0 , composition)
    Cv = PropsSI('CVMASS' , 'T' , T0 , 'P' , p0 , composition)
    n = Cp / Cv  # показатель адиабаты = показатель политропы 1.25 в некоторых источниках
    Tk = T0 * (1 + 1 / kpd_ad * ((nk) ** ((n - 1) / n) - 1))
    Tk0 = 0
    while abs(Tk0 - Tk) > 0.1:
        Cp = PropsSI('C' , 'T' , (T0 + Tk) / 2 , 'P' , p0 , composition)
        Cv = PropsSI('CVMASS' , 'T' , (T0 + Tk) / 2 , 'P' , p0 , composition)
        n = Cp / Cv
        Tk0 = Tk
        Tk = T0 * (1 + 1 / kpd_ad * ((nk) ** ((n - 1) / n) - 1))

    lk = 0
    for i in range(n_kol_stup_komp):
        lki = n / (n - 1) * R * T0 * ((nk) ** ((n - 1) / n) - 1) / kpd_ad  # политропное сжатие - рабочий процесс
        lk = lk + lki
    Qk = Cp * (Tk - T0)*n_kol_stup_komp
    print("Температура в конце сжатия " , round(Tk , 2) , "К")
    print("Удельная работа компрессора " , round(lk / 1000 , 2) , "кДж/кг")
    print("Теплота сжатия " , round(Qk / 1000 , 2) , "кДж/кг")
    return lk , Tk , Qk

def detander(p_in , T_in , p_out): #T_out , Ld
    k_det = 0.75
    s_in = PropsSI('S' , 'T' , T_in , 'P' , p_in , 'REFPROP::Air')
    h_in = PropsSI('H' , 'S' , s_in , 'P' , p_in , 'REFPROP::Air')
    h_out_t = PropsSI('H' , 'S' , s_in , 'P' , p_out , 'REFPROP::Air') #    s_out = s_in
    Ld = k_det * (h_in - h_out_t)
    T_out = PropsSI('T' , 'H' , h_in - Ld , 'P' , p_out , 'REFPROP::Air')
    print("Температура на выходе из детандера = " , round(T_out , 2) , " К")
    print("Удельная работа детандера = " , round(Ld / 1000 , 2) , " кДж/кг")
    return T_out , Ld

def detanderORC(p_in , T_in , p_out,Toc): #T_out , Ld
    k_det = 0.75
    s_in = PropsSI('S' , 'T' , T_in , 'P' , p_in , 'R134')
    h_in = PropsSI('H' , 'S' , s_in , 'P' , p_in , 'R134')
    s_out = s_in
    h_out_t = PropsSI('H' , 'S' , s_out , 'P' , p_out , 'R134')
    Ld = k_det * (h_in - h_out_t)
    h_out = h_in - Ld
    T_out = PropsSI('T' , 'H' , h_out , 'P' , p_out , 'R134')
    s_out = PropsSI('S' , 'T' , T_out , 'P' , p_out , 'R134')
    exergy_td = ((h_in - h_out) - Toc * (s_in - s_out))
    lossExergy_td = exergy_td - Ld
    E_td = 1 - lossExergy_td / Ld
    print("Потери эксергии в детандере" , lossExergy_td)
    print("Эксергетический КПД детандера" , E_td)
    print("Температура на выходе из детандера = " , round(T_out , 2) , " К")
    print("Удельная работа детандера = " , round(Ld / 1000 , 2) , " кДж/кг")
    return T_out , Ld

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

def Pump(p0 , T0 , pk, composition ): #Lp , Tp
    k_pump = 0.75
    density = PropsSI("D" , "T" , T0 , "P|liquid" , (p0 + pk) / 2 , composition)
    Lp = (pk - p0) / density / k_pump
    h0 = PropsSI("H" , "T" , T0 , "P|liquid" , p0 , composition)
    s0 = PropsSI("S" , "T" , T0 , "P|liquid" , p0 , composition)
    hp = h0 + Lp
    Tp = PropsSI("T" , "H" , hp , "P" , pk , composition)
    sp = PropsSI("S" , "T" , Tp , "P|liquid" , pk , composition)
    Exergy_p = ((hp - h0) - T0 * (sp - s0))
    lossExergy_p = Lp - Exergy_p
    E_p = 1 - lossExergy_p / Lp
    print("--Эксергетический КПД насоса" , E_p)
    # print("Температура после насоса " , round(Tp , 2) , "К")
    # print("Работа насоса " , round(Lp / 1000 , 2) , " кДж/кг")
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

def TurbinaORC(Toc , Qkm , Gtur , pt1 , Tt1 , pt2 , Tkm, st ): #Ld_sum*k_tur, Qkm
    n_st_det = 1
    while (pt1 / pt2) ** (1 / n_st_det) > st:
        n_st_det += 1

    k_tur = 0.75
    nd = (pt1 / pt2) ** (1 / n_st_det)
    print("Степень расширения " , nd , " Количество ступеней расширения " , n_st_det)
    Ld_sum = 0
    for i in arange(n_st_det):
        pt2 = pt1 / nd
        Cp = PropsSI('CPMASS' , 'T' , Tt1 , 'P' , (pt1+pt2)/2 , 'R134')
        T_heat = Qkm/Cp/Gtur+Tt1
        Theat0 = 0
        while (T_heat - Theat0)>0.1:
            Theat0 = T_heat
            Cp = PropsSI('CPMASS' , 'T' , (Tt1+T_heat)/2 , 'P' , (pt1+pt2)/2 , 'R134')
            T_heat = Qkm / Cp / Gtur + Tt1
        if T_heat < Toc:
            T_heat = Toc
        if T_heat > Tkm:
            T_heat = Tkm - 5
            if T_heat > 550:  # ограничение температуры по термомаслу
                T_heat = 550
        print("Температура подогретого воздуха " , T_heat , " К")
        Qkm = Qkm - Cp * Gtur*(T_heat-Tt1)
        print("Давление стало: " , pt2 / 1000000 , "МПа")
        Td , Ld = detander(pt1 , T_heat , pt2,Toc)
        Ld_sum = Ld_sum + Ld
        Tt1 = Td
        pt1 = pt2
    return Ld_sum*k_tur, Qkm

def maxTempbeforeTrottle(p_max , p_trot):
    H_liq = PropsSI('H' , 'P' , p_trot , 'Q' , 1 , 'REFPROP::Air')
    H_max = H_liq
    T_max = PropsSI('T' , 'P' , p_max , 'H' , H_max , 'REFPROP::Air')
    return T_max

def ChartTtoL():
    masLd1 = []
    masLd2 = []
    masLd3 = []
    masT = []
    pd1 = 6000000
    pd2 = 10000000
    pd3 = 18000000
    p0 = 101325
    for T in arange(273 , 773 , 1):
        Td1 , Ld1 = detander(pd1 , T , p0,300)
        Td2 , Ld2 = detander(pd2 , T , p0,300)
        Td3 , Ld3 = detander(pd3 , T , p0,300)
        masT.append(T)
        masLd1.append(Ld1 / 1000)
        masLd2.append(Ld2 / 1000)
        masLd3.append(Ld3 / 1000)

    x = np.array(masT)
    y1 = np.array(masLd1)
    y2 = np.array(masLd2)
    y3 = np.array(masLd3)
    plt.title("")
    plt.xlabel("Температура на входе в турбину, К")
    plt.ylabel("Удельная работа турбины, кДж/кг")
    plt.plot(x , y1 , color="green" , label="pk = 6 МПа")
    plt.plot(x , y2 , color="red" , label="pk = 10 МПа")
    plt.plot(x , y3 , color="blue", label="pk = 18 МПа")
    plt.legend()
    plt.grid(True)
    plt.show()

def ChartPtoL():
    masLd1 = []
    masLd2 = []
    masLd3 = []
    masPd = []
    T1 = 300
    T2 = 400
    T3 = 500
    p0 = 101325
    for pd in arange(1*10**6 , 20*10**6 , 500000):
        Td1 , Ld1 = detander(pd , T1 , p0,300)
        Td2 , Ld2 = detander(pd , T2 , p0,300)
        Td3 , Ld3 = detander(pd , T3 , p0,300)
        masPd.append(pd/1000000)
        masLd1.append(Ld1 / 1000)
        masLd2.append(Ld2 / 1000)
        masLd3.append(Ld3 / 1000)

    x = np.array(masPd)
    y1 = np.array(masLd1)
    y2 = np.array(masLd2)
    y3 = np.array(masLd3)
    plt.title("")
    plt.xlabel("Давление на входе в турбину, МПа")
    plt.ylabel("Удельная работа турбины, кДж/кг")
    plt.plot(x , y1 , color="green" , label="Т = 300 К")
    plt.plot(x , y2 , color="red" , label="Т = 400 К")
    plt.plot(x , y3 , color="blue", label="Т = 500 К")
    plt.legend()
    plt.grid(True)
    plt.show()

