#Скрипт для полностью автоматизированного расчёта pKa в ORCA 5.
Рассчёт выполняется через вычисление dG реакции кислотно-основного равновесия исследуемого соединения с эталоном.

AH+REF-=A-+REFH --->dG  
`pKa=pKa(REF)+LOG10(e)*dG*2.6e6/(R*T)`

Для вычисления dG используется комбинированый подход:

dG=dG(r2SCAN-3c)-SP(r2SCAN-3c)+SP(revDSD-PBEP86-D4/aug-cc-pVTZ).

В качестве REF используется DMSO (pKa 35), раствор моделируется CPCM(DMSO)

Перед использованием нужно проверить параметры указанные в скрипте:

ORCA=os.path.expanduser('~')+"/orca/orca" #путь к ORCA  
SOLVENT="DMSO" #Используемый растворитель  
MAXCORE="12000" # Обём памяти на ядро.  
NPROC="32" # Количество используемых ядер.  

##Использование:

`$automatic_calc_pka.py AH.xyz A-.xyz [BH.xyz B-.xyz]`

Все промежуточно вычесленные энергии пишутся в файл r2scan-3c_opt_freq_revDSDPBEP86d4_sp_report.csv, вычесленные значения pKa пишутся в файл automatic_pka_report_file.csv

##Зависимости:

>- ORCA 5.0.4  
>- python >3.8  
>- numpy  

__________

Александр Белый 
2023
ИОХ РАН 
Лаб. №6 Химии диазосоединений.
