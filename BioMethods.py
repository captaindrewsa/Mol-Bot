import re
import math
import requests
from config import TOKEN
from random import randint

#Константы для соседей
NEAREST_NEIGHBORS= {'AA':[[-9.1,-0.0240],[-6.6,-0.0184]],
'AT':[[-8.6,-0.0239],[-5.7,-0.0155]],
'TA':[[-6.0,-0.0169],[-8.1,-0.0226]],
'CA':[[-5.8,-0.0129],[-10.5,-0.0278]],
'GT':[[-6.5,-0.0173],[-10.2,-0.0262]],
'CT':[[-7.8,-0.0208],[-7.6,-0.0192]],
'GA':[[-5.6,-0.0135],[-13.3,-0.0355]],
'CG':[[-11.9,-0.0278],[-8.0,-0.0194]],
'GC':[[-11.1,-0.0267],[-14.2,-0.0349]],
'GG':[[-11.0,-0.0266],[-12.2,-0.0297]],
'TT':[[-9.1,-0.0240],[-6.6,-0.0184]],
'TG':[[-5.8,-0.0129],[-10.5,-0.0278]],
'AC':[[-6.5,-0.0173],[-10.2,-0.0262]],
'AG':[[-7.8,-0.0208],[-7.6,-0.0192]],
'TC':[[-5.6,-0.0135],[-13.3,-0.0355]],
'CC':[[-11.0,-0.0266],[-12.2,-0.0297]],
'AU':[[-8.6,-0.0239],[-5.7,-0.0155]], #дальше дополнительное для РНК
'UA':[[-6.0,-0.0169],[-8.1,-0.0226]],
'GU':[[-6.5,-0.0173],[-10.2,-0.0262]],
'CU':[[-7.8,-0.0208],[-7.6,-0.0192]],
'UU':[[-9.1,-0.0240],[-6.6,-0.0184]],
'UG':[[-5.8,-0.0129],[-10.5,-0.0278]],
'UC':[[-5.6,-0.0135],[-13.3,-0.0355]],}


'''
Константа (учитывает зарождение спирали. во время отжига / плавления) 
A = -0.0108
Газовая постоянная
R = 0.00199
Концентрация олигонуклеотида в M или моль л-1 
(Я использую 0,0000005, т.е. 0,5 мкМ) 
C = 0.0000005
Концентрация ионов натрия в М или моль л -1 
(Я использую 0.05, т.е. 50 мМ)
CNa = 0.05
'''
#Вспомогательные функции
def CountDeltaHandDeltaSDNA(NA):
    deltaH=0
    deltaS=0
    for idx in range(0,len(NA)-1):
        deltaH += NEAREST_NEIGHBORS[NA[idx]+NA[idx+1]][0][0]
        print(NEAREST_NEIGHBORS[NA[idx]+NA[idx+1]][0][0])
        deltaS += NEAREST_NEIGHBORS[NA[idx]+NA[idx+1]][0][1]
    return [deltaH,deltaS]

def CountDeltaHandDeltaSRNA(NA):
    deltaH=0
    deltaS=0
    for idx in range(0,len(NA)-1):
        deltaH += NEAREST_NEIGHBORS[NA[idx]+NA[idx+1]][1][0]
        deltaS += NEAREST_NEIGHBORS[NA[idx]+NA[idx+1]][1][1]
    return [deltaH,deltaS]

def BeautifulOutput(DISCorDATA):
    out_string = '<b>'
    i = 1
    for elem in DISCorDATA:
        out_string+='{}) {}\n'.format(i,elem)
        i+=1
    out_string+='</b>'
    return out_string


#Разворот последовательности
def ReverseNa(NA):
    ret_NA=''
    
    if len(NA)%70 == 0:
        for i in range(0,len(NA)//70):
            ret_NA+=(NA[::-1])[70*i:70*(i+1)]
            ret_NA+='\n'
        return ret_NA.upper()
    else:
        for i in range(0,len(NA)//70):
            ret_NA+=(NA[::-1])[70*i:70*(i+1)]
            ret_NA+='\n'
        ret_NA+=NA[::-1][len(NA)-len(NA)%70:]+'\n'
        return ret_NA.upper()

    



#Вывод комплиментарной цепочки
def ComplDna(NA):
    COMPL_DNA = ''
    for i in NA.upper():
        if i == 'A':
            COMPL_DNA += 'T'
        elif (i == 'T') | (i == "U"):
            COMPL_DNA += 'A'
        elif i == 'G':
            COMPL_DNA += 'C'
        elif i == 'C':
            COMPL_DNA += 'G'
        elif i == '\n':
            COMPL_DNA += '\n'
        else:
            return 'Ошибка в нуклеотиде {}'.format(len(COMPL_DNA))
        if len(COMPL_DNA)%70 == 0:
            COMPL_DNA+='\n'
    return COMPL_DNA

#Если надо в РНК
def ComplRna(NA):
    COMPL_RNA = ''
    for i in NA.upper():
        if i == 'A':
            COMPL_RNA += 'U'
        elif (i == 'T') | (i == 'U'):
            COMPL_RNA += 'A'
        elif i == 'G':
            COMPL_RNA += 'C'
        elif i == 'C':
            COMPL_RNA += 'G'
        else:
            return 'Ошибка в нуклеотиде {}'.format(len(COMPL_RNA))
        if len(COMPL_RNA)%70 == 0:
            COMPL_RNA+='\n'
    return COMPL_RNA


#Темепратура отжига олигов
def TemOlig(NA):
    if 15<=len(NA)<=120:
        TOligDNA = (CountDeltaHandDeltaSDNA(NA)[0]/(-0.0108 + CountDeltaHandDeltaSDNA(NA)[1] + (0.00199*math.log(0.0000005/4.0)))) - 273.15+16.6*math.log10(0.05)
        TOligRNA = (CountDeltaHandDeltaSRNA(NA)[0]/(-0.0108 + CountDeltaHandDeltaSRNA(NA)[1] + (0.00199*math.log(0.0000005/4.0)))) - 273.15+16.6*math.log10(0.05)
        return [TOligDNA,TOligRNA]
    else:
        TOligDNA= 2*(len(re.findall("A",NA))+len(re.findall("T",NA)))+4*(len(re.findall("C",NA))+len(re.findall('G',NA)))-7
        TOligRNA= 2*(len(re.findall("A",NA))+len(re.findall("U",NA)))+4*(len(re.findall("C",NA))+len(re.findall('G',NA)))-7
        return [TOligDNA,TOligRNA] 

#Зааворачивание в фасту
def OutFasta(file_name,description, NA):
    with open('fasta_files/{}.fasta'.format(file_name), 'w+', encoding='utf-8') as fasta:
        NA_lost=NA
        fasta.write('{}\n'.format(description))
        while True:
            if len(NA_lost)==0:
                break
            else:
                fasta.write(NA_lost[:70]+'\n')
                NA_lost=NA_lost[70:]
                continue

#Принятие и обработка файла
    #ДОДЕЛАТЬ определяет фаста или мультифаста
    #Например РАЗВЕРНУТЬ ВСЕ, КОМПЛИМЕНТАРНОЕ ВСЕ, ТЕМПЕРАТУРУ ВСЕГО, ПРАЙМЕРЫ И ЕЩЕ ЧЕНИТЬ (но это уже после tkinter, наверное)

def ListDataInFasta(file_info):
    with open('TEST.txt', 'r', encoding='utf-8') as f:
        lines = f.readlines()
        DISCRIPTIONS =  lines[::2]
        DATA =  lines[1::2]
        return [DISCRIPTIONS, DATA]

#Lambda-red recombination
def LambdaRedRecombination(mRNA, genomeInsertion, start, stop):
    #Подбор олигонуклеотидов для рекомбинации
    gDNA_53 = ComplDna(mRNA.strip())

    gDNA_olig1_35 = ComplDna(gDNA_53[:20])
    gDNA_olig2_53 = ((gDNA_53[::-1])[:20])[::-1]

    genome35 = ComplDna(genomeInsertion)

    genome_olig1_35 = (genome35[start-51:start-1])
    genome_olig2_53 = (genomeInsertion[stop+50:stop:-1])[::-1]

    homology_olig1_53 = (genome_olig1_35+gDNA_olig1_35)[::-1] #Ретёрнить
    homology_olig2_53 = gDNA_olig2_53+genome_olig2_53         #Ретёрнить
    #Подбор двух пар праймеров для проверки
    #Праймеры на вставку
    check_insert_primer1 = searchOligsInsertion(genomeInsertion, start, stop)[0]
    check_insert_primer2 = searchOligsInsertion(genomeInsertion, start, stop)[1]
    #праймеры на wrong way
    check_wrongway_primer1 = check_insert_primer1
    check_wrongway_primer2 = ((gDNA_53[::-1])[:len(check_wrongway_primer1)])[::-1] #Ретёрнуть

    return [[homology_olig1_53, homology_olig2_53],[check_insert_primer1, check_insert_primer2],[check_wrongway_primer1,check_wrongway_primer2]]



#Поиск олигонуклеотидов на вставку
def searchOligsInsertion(genomeInsertion, start, stop):
    zone53 = genomeInsertion[start-100:start-25]
    zone35 = ComplDna(genomeInsertion)[stop+25:stop+100]

    for i in range(0,75):
        temp1 = TemOlig(zone53[i:i+26])[0]
        for e in range(0,75):
            temp2 = TemOlig(zone35[e:e+26])[0]
            if math.fabs((temp1-temp2))<=5:
                return [zone53[i:i+26], zone35[e:e+26][::-1]]
            else:
                continue
#Поиск олигонуклеотидов на WrongWay
def searchOligsWrongWay(mRNA, genomeInsertion):
    pass

#Вставка гена по гомологии
def HomologyRecombination(GeneWithShoulder, genomeInsertion):
    Shoulder1=''
    Shoulder2=''
    i=1
    #Ищем плечо 5-3:
    while True:
        if genomeInsertion.find(GeneWithShoulder[:i]) != (-1):
            i+=1
            continue
        else:
            Shoulder1 = GeneWithShoulder[:i]
            i = 1
            break
    #Ищем плечо 3-5:
    while True:
        if (genomeInsertion[::-1]).find(GeneWithShoulder[::-1][:i]) != (-1):
            i+=1
            continue
        else:
            Shoulder2 = (GeneWithShoulder[::-1][:i])
            i = 1
            break
    #Вставляем в геном наш ген по плечам
    idx_start = genomeInsertion.find(Shoulder1)
    idx_end = (genomeInsertion.find(Shoulder2)) + len(Shoulder2)
    final_genome ='{0}{1}{2}'.format(genomeInsertion[:idx_start],GeneWithShoulder,genomeInsertion[idx_end:])
    return final_genome









#Выдача карты рестрикции с номерами сайтов
#Обрезание по сайтам рестрикции с выдачей фаст.
#Создание фаст нужных олигов с нужными липкими концами
#Создание фаст для необходимой последовательности с плечами гомологии в нужном месте
#ФОРЕЗ???
#Подобрать праймеры для "colony PCR"

