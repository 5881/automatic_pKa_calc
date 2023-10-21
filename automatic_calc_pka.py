#!/bin/python
#Александр Белый 
#2023
#ИОХ РАН 
#Лаб. №6 Химии диазосоединений.
#Alexander Belyy IOC RAS
#Скрипт автоматического рассчёта pKa исходя из структуры кислоты и аниона
#Для нахождения pKa используется линейная связь dG реакции исследуемой 
#кислоты со стандартом (DMSO)
#на первом этапе производится оптимизация структур методом r2scan-3c
#и вычисление частот, восле чего энергия SP уточняется двойным гибридом
#revDSD-PBEP86(d4)/aug-cc-PVTZ
#Полученые гибридные значения dG используются для рассчётов
#dG для DMSO рассчитана заранее и берётся в виде константы
#usage:
#$automatic_calc_pka.py acid1 anion1 acid2 anion2 ...
#
from collections import namedtuple
import os,sys
import numpy as np
ORCA=os.path.expanduser('~')+"/orca/orca"
ENERGYREPORTFILE="r2scan-3c_opt_freq_revDSDPBEP86d4_sp_report.csv"
PKAREPORTFILE="automatic_pka_report_file.csv"
SOLVENT="DMSO"
MAXCORE="12000"
NPROC="32"
CONTINUE=True

DGREF=0.5241152317099704

T=298
R=8.314
PKAREF=35 #pKa стандарта
KAREF=pow(10,-PKAREF)
LOG10EXP=0.43429448190325182765

def calc_pka_from_dg(dg):
	return PKAREF+LOG10EXP*dg*2.6e6/(R*T)

Energy = namedtuple('Energy' , 'sp g h gcor hcor geomok')
#Combinedenergy = namedtuple('Combinedtnergy' , 'sp g h gcor hcor geomok')

def copyfile(src, dst):
	cmd="cp -R "+src+" "+dst
	os.system(cmd)

def getxyz_from_file(filename):
	xyz=""
	with open(filename, 'r', encoding="utf-8", errors="ignore") as xyz_file:
		for line in xyz_file.readlines():
			if (len(line.split()) == 4 and line.split()[1][-1].isdigit()):
				xyz+=line
	return xyz

def output_terminate_status(outputfile):
	with open(outputfile, 'r') as out_file:
		for s in out_file:
			if '****ORCA TERMINATED NORMALLY****' in s:
				return True
	return False

def output_freq_analyse(outputfile):
	#Теперь анализируем всю термохимию в out файле, так удобнее
	gibbs_done=0 
	vib_freq_done=0
	geom_ok=1
	G=0
	H=0
	G_el=0
	E=0
	Ecor=0
	with open(outputfile, 'r') as out_file:
		for s in out_file:
			if 'FINAL SINGLE POINT ENERGY' in s: E=float(s.split()[-1])
			if 'VIBRATIONAL FREQUENCIES' in s: vib_freq_done=1
			if vib_freq_done and '***imaginary mode***' in s:
				print("Мнимые частоты, геометрия возможно неверна")
				geom_ok=0
			if 'Final Gibbs free energy' in s: G=float(s.split()[-2])
			if 'G-E(el)' in s: G_el=float(s.split()[-4])
			if 'Total correction' in s: Ecor=float(s.split()[2])
	print(outputfile, "E=", E, "G=", G,"Ecor=",Ecor,"G-el=",G_el, "OK=",geom_ok)
	#print("eee",E+G_el)
	#return E, G, Ecor, G_el, geom_ok
	#Energy = namedtuple('Energy' , 'sp g h gcor hcor geomok')
	H=E+Ecor
	return Energy(E,G,H,G_el,Ecor,geom_ok)

def method_run_xyzfile(job,method_str,config,xyzfile):
	#NPROC="12"
	"""
	заготовка метода, чтоб не плодить сущьности
	job - добавляется в качестве суффикса в имя файла
	method_str - задаёт метод базисс, модель растворителя и тп
	#!pbeh-3c freq TightSCF TightOPT cpcm(SOLVENT) xyzfile
	config - дополнительные параметры
	#%geom MaxIter 500 end
	xyzfile - файл с координатами
	"""
	OPTJOB=False
	RUN=True
	if "opt" in method_str.lower():OPTJOB=True
	
	#job="_pbeh3c_opt_freq"
	basename=xyzfile.replace('.xyz', '').rsplit('_')[0]
	basename=basename+job
	dirname=basename
	inpfile=basename+".inp"
	outputfile=basename+".out"
	xyzoutputfile=basename+".xyz"
	charge=0
	xyz=getxyz_from_file(xyzfile)
	if(basename.find("-")!=-1): charge=-1
	if(basename.find("+")!=-1): charge=1
	if(os.path.isdir(dirname) is False):
		os.mkdir(dirname)
	os.chdir(dirname)
	
	if(os.path.isfile(outputfile) is True):
		if(output_terminate_status(outputfile) is False):
			if(CONTINUE and OPTJOB): xyz=getxyz_from_file(xyzoutputfile)
			#cmd='find . -not -name "'+xyzoutputfile+'" -delete'
			os.system("rm *")
		else: RUN=False
	
	method_str=method_str.replace('SOLVENT',SOLVENT)
	if(RUN is True):
		with open(inpfile, 'w') as inp_file:
			print('%pal nprocs {0} end'.format(NPROC),file=inp_file)
			print('%maxcore {0}'.format(MAXCORE),file=inp_file)
			print(method_str,file=inp_file)
			print(config,file=inp_file)
			print('* xyz {0} 1'.format(charge),file=inp_file)
			print(xyz,"*",sep='',file=inp_file)
		orcacmd=ORCA+" "+inpfile+"|tee "+outputfile
		os.system(orcacmd)
	#финальную геометрию лучше сохранить, если это оптимизация
	if OPTJOB: copyfile(xyzoutputfile,"..")
	#нововведение, используем именованый кортеж для хранения энергий
	energys=output_freq_analyse(outputfile)
	#outputfile=dirname+"/"+outputfile
	#xyzfile=dirname+"/"+xyzoutputfile
	os.chdir("..")
	if OPTJOB: return energys,xyzoutputfile
	return energys

def r2scan_opt_freq(xyzfile):
	job="_r2scan_opt_freq"
	method_str="!r2scan-3c freq TightSCF veryTightOPT defgrid3 cpcm(SOLVENT) xyzfile"
	config="%geom MaxIter 500 end"
	return method_run_xyzfile(job,method_str,config,xyzfile)

def revdsd_pbep86_d4(xyzfile):
	#global NPROC="12"
	job="_revDSD_PBEP86_D4"
	method_str="!D4 aug-cc-pvtz rijk aug-cc-pvtz/C aug-cc-pvtz/jk TightSCF cpcm(SOLVENT)"
	config="""
%basis
PCDTrimBas Overlap # Trim the orbital basis in the overlap metric
#PCDTrimAuxJ Coulomb # Trim the AuxJ basis in the Coulomb metric
PCDTrimAuxJK Coulomb # Trim the AuxJK basis in the Coulomb metric
PCDTrimAuxC Coulomb # Trim the AuxC basis in the Coulomb metric
PCDThresh -1 # Threshold for the PCD: chosen automatically if <0
end
%method
   Exchange X_PBE
   Correlation C_P86
   ScalHFX 0.69
   ScalDFX 0.31
   ScalGGAC 0.4210
   ScalLDAC 0.4210
   ScalMP2C 1.0
   LDAOpt C_VWN5
   D3S6 0.5132
   D3S8 0.0
   D3A1 0.44
   D3A2 3.60
end
%mp2
   DoSCS true
   PS 0.5922
   PT 0.0636
end
"""
	return method_run_xyzfile(job,method_str,config,xyzfile)

def calc_combined_g(xyzfile):
	r2scan_energy,xyzoutputfile=r2scan_opt_freq(xyzfile)
	dh_sp_energy=revdsd_pbep86_d4(xyzoutputfile)
	G=r2scan_energy.gcor+dh_sp_energy.sp
	H=r2scan_energy.hcor+dh_sp_energy.sp
	#главное застолбить промежуточный результат
	if(os.path.isfile(ENERGYREPORTFILE) is False):
		with open(ENERGYREPORTFILE, 'w') as reportfile:
			print(
			"xyzfile, geomok, r2scan-3c(sp), r2scan-3c(H), r2scan-3c(G), hcor, gcor, revDSD-PBEP86(d4)(sp), H, G"
			,file=reportfile)
			#print("filename,E,G,H,Hcor,Gcor,geomok",file=reportfile)
	with open(ENERGYREPORTFILE, 'a') as reportfile:
		print('{0}, {1}, {2}, {3}, {4}, {5}, {6}, {7}, {8}, {9}'.format(
			xyzfile,
			r2scan_energy.geomok,
			r2scan_energy.sp,
			r2scan_energy.h,
			r2scan_energy.g,
			r2scan_energy.hcor,
			r2scan_energy.gcor,
			dh_sp_energy.sp,
			H,
			G),file=reportfile)
	return Energy(dh_sp_energy.sp,G,H,r2scan_energy.gcor,r2scan_energy.hcor,r2scan_energy.geomok)

def calc_dg(reactantxyz,productxyz):
	energy_reactant=calc_combined_g(reactantxyz)
	energy_product=calc_combined_g(productxyz)
	dG=energy_product.g-energy_reactant.g
	return dG

def calc_pka(reactantxyz,productxyz):
	energy_reactant=calc_combined_g(reactantxyz)
	energy_product=calc_combined_g(productxyz)
	dG_pair=energy_product.g-energy_reactant.g
	dG_reaction=dG_pair-DGREF
	pKa=calc_pka_from_dg(dG_reaction)
	#Обязательно всё сохранить
	print("dG(AH>A-)={0}\npKa(HA)={1}".format(dG_reaction,pKa))
	if(os.path.isfile(PKAREPORTFILE) is False):
		with open(PKAREPORTFILE, 'w') as reportfile:
			print(
			"acid xyzfile, anion xyzfile, dG, dG-dG(REF),pKa"
			,file=reportfile)
			#print("filename,E,G,H,Hcor,Gcor,geomok",file=reportfile)
	with open(PKAREPORTFILE, 'a') as reportfile:
		print('{0}, {1}, {2}, {3}, {4}'.format(
			reactantxyz,
			productxyz,
			dG_pair,
			dG_reaction,
			pKa),file=reportfile)


N=int(len(sys.argv[1:])/2)
for i in range(N):
	calc_pka(sys.argv[2*i+1],sys.argv[2*i+2])

#print("dGREF=",DGREF)
