#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
""""""""""""""""""""""
PythonMeta:
Meta-Analysis Package
By dhy @ SHUTCM 2017
"""""""""""""""""""""""
import os
import math

import matplotlib
import matplotlib.pyplot as plt
from pylab import gca

import logging
from os import path

logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s %(name)s %(levelname)s %(message)s',
                    filename=path.join(path.normpath(path.dirname(__file__)), '','meta.log'),
                    datefmt='[%d/%b/%Y %H:%M:%S]',
                    filemode='a')
logger = logging.getLogger(__name__)
import scipy.stats

_INFO_ = """
PythonMeta
=====================

Name = PythonMeta
Version = 1.11 
Author = Deng Hongyong (dhy)
Email = dephew@126.com
URL = www.pymeta.com

This is a Meta-Analysis package. 

This module was designed to perform some Evidence-based medicine (EBM) tasks, such as:

* Combining effect measures (OR, RR, RD for count data and MD, SMD for continuous data);
* Heterogeneity test(Q/Chi-square test);
* Subgroup analysis;
* Plots drawing: forest plot, funnel plot, etc.

Statistical algorithms in this software cited from:
**Jonathan J Deeks and Julian PT Higgins, on behalf of the Statistical Methods Group of The Cochrane Collaboration. Statistical algorithms in Review Manager 5, August 2010.**

Please cite me in any publictions like:
**Deng Hongyong. PyMeta, Python module of Meta-analysis, cited 20xx-xx-xx (or your time); 1 screen(s). Available from URL: http://www.pymeta.com**

This is an ongoing project, so, any questions and suggestions from you are very welcome.

Please visit www.pymeta.com for samples and more details of usage.
"""

def Help():
    print(_INFO_)

#class of data processing
class Data:
    def __init__(self):
        self.datatype=''    #data type:'CATE','CONT'
        self.studies=[]     #all studies
        self.subgroup=[0]   #[0]number of subgroups,[1]... name of subgroups,study id...
        self.nototal=False  #flag of do NOT calculate the total effect
        logger.info("Load Data().")
        
    #get studies from lines
    #input:formatted lines, e.g. studies.txt
    #output: studies[[],...]
    def getdata(self,lines):
        self.datatype=chkdtype(self.datatype)

        def _is_number(s):
            try:
                float(s)
                return True
            except :
                return False

        sty=[];rult=[];subgrp=[];s=""
        if self.datatype=='CATE':
            stylen=5
            def _getnum(s):
                return int(s)
        elif self.datatype=='CONT':
            stylen=7
            def _getnum(s):
                return float(s)
        else:
            Err="error.(no datatype setting:'CATE', 'CONT')"
            raise Exception(Err)

        for line in lines:
            line=line.rstrip("\n").rstrip("\r").strip() #trim
            line=line.replace("；",",").replace(";",",")
            line=line.replace("，", ",").replace("\t", ",")
            line=line.replace("\'", "").replace("\"", "").replace(" ","")
            line=line.replace(",,", ",").replace(",,",",")
            line=line.strip(",").strip()
            if line[0:1]=='#':
                continue
            elif len(line)>9 and line.count(",")==stylen-1 :
                sty=line.split(",")
                j=0;sty0=[0]*(stylen-1);sty0.append('')
                for i in range(stylen):
                    if _is_number(sty[i]) :
                        sty0[j]=_getnum(sty[i])
                        j+=1
                    elif len(sty[i])>0 :
                        sty0[stylen-1]=sty[i].strip()[0:20]  #study name

                if j==stylen-1 and len(sty0[stylen-1])>0 :
                    rult.append(sty0)
                    if len(rult)>200:
                        Err="error.(too many studies(<200))"
                        raise Exception(Err)
                    subgrp.append(len(rult)-1)
            elif line[0:12].lower()=='<comparison>': #not use here
                self.comparison=midstr(line,"=")
            elif line[0:9].lower()=='<outcome>':
                self.outcome=midstr(line,"=")
            elif line[0:9].lower()=='<nototal>':
                self.nototal=True
            elif line[0:10].lower()=='<subgroup>':
                self.subgroup[0]+=1
                if self.subgroup[0]>20:
                    Err="error.(too many subgroups(<20))"
                    raise Exception(Err)
                self.subgroup.append([midstr(line,"=") if len(midstr(line,"="))>0 else 'subgroup_{0:d}'.format(self.subgroup[0])]+subgrp)
                subgrp=[]

        if len(rult)<1:
            Err="error.(no data loaded)"
            raise Exception(Err)

        self.studies=rult
        return rult

    #load text file to studies
    #input:filename like "c:\\1.txt"
    #output: lines
    def readfile(self,fname='studies.txt'):
        try:
            f = open(fname,"r")
            lines = f.readlines()   
        except :
            Err="error.(no data file or file name error)"
            raise Exception(Err)

        f.close()
        if len(lines)<1:
            Err="error.(no data loaded from file)"
            raise Exception(Err)
        return lines

#class of Meta-analysis
class Meta:
    def __init__(self):
        self.datatype=''       #'CATE': binary/categorical/dichotomous data; and 'CONT': continuous data
        self.models='Fixed'    #'Fixed';'Random'
        self.effect=''         #'OR':odds ratio; 'RR': risk ratio; 'RD':risk difference; 'MD':weighted mean diff; 'SMD':standard mean diff
        self.algorithm=''      #'MH':Mantel-Haenszel;'Peto';'IV':Inverse variance;'IV-Heg'(DEFAULT),'IV-Cnd','IV-Gls':for SMD algorithms
        self.studies=[]        #all studies
        self.subgroup=[0]      #[0]number of subgroups,[1]:[name of subgroup1, study id,...],[n]...
        self.nototal=False     #force to do NOT calculate overall effect
        logger.info("Load Meta().")
        
    #check methods
    def chkmethods(self):
        self.datatype=chkdtype(self.datatype)
        md=self.models.lower().replace("-"," ").replace("_"," ").strip()
        md=md.split(" ")[0]
        a_fix=["fix","fixed"]
        a_rnd=["rnd","rand","random"]
        if md in a_fix:
            self.models="Fixed"
        elif md in a_rnd:
            self.models="Random"
        else:
            Err="error.(effect models should be:'Fixed' or 'Random')"
            raise Exception(Err)

        md=self.effect.lower().replace("-"," ").replace("_"," ").strip()
        a_or=["or","odds","odds ratio"]
        a_rr=["rr","risk ratio"]
        a_rd=["rd","risk difference","risk diff"]
        a_md=["md","wmd","mean difference","weighted mean difference","mean diff","weighted mean diff"]
        a_smd=["smd","standard mean difference","standard mean diff"]
        if md in a_or:
            self.effect="OR"
        elif md in a_rr:
            self.effect="RR"
        elif md in a_rd:
            self.effect="RD"
        elif md in a_md:
            self.effect="MD"
        elif md in a_smd:
            self.effect="SMD"
        else:
            Err="error.(effect measures should be:'OR', 'RR', 'RD', 'MD' or 'SMD')"
            raise Exception(Err)

        if self.datatype=='CATE':
            if not self.effect in "RR,OR,RD":
                Err="error.(bad eff_size in CATE)"
                raise Exception(Err)
        elif self.datatype=='CONT':
            if not self.effect in "MD,SMD":
                Err="error.(bad eff_size in CONT)"
                raise Exception(Err)

        md=self.algorithm.lower().replace("-"," ").replace("_"," ").replace(","," ").strip()
        a_mh=["mh", "m h", "mantel haenszel"]
        a_pt=["peto"]
        a_iv=["iv","i v","inverse variance"]
        a_ivh=["iv heg","heg","hedges","hedges' adjusted g"]
        a_ivc=["iv cnd","cnd","cohen","cohens","cohen's","cohen's d"]
        a_ivg=["iv gls","gls","glass","glass's","glass's d"]
        if md in a_mh:
            self.algorithm="MH"
        elif md in a_pt:
            self.algorithm="Peto"
        elif md in a_iv:
            self.algorithm="IV"
        elif md in a_ivh:
            self.algorithm="IV-Heg"
        elif md in a_ivc:
            self.algorithm="IV-Cnd"
        elif md in a_ivg:
            self.algorithm="IV-Gls"
        else:
            Err="error.(algorithm should be:'MH', 'Peto', 'IV', 'IV-Heg', 'IV-Cnd' or 'IV-Gls')"
            raise Exception(Err)

        if self.effect=='OR':
            if not self.algorithm in "MH,Peto,IV":
                Err="error.(algorithm for OR should be:'MH', 'Peto' or 'IV')"
                raise Exception(Err)
        elif self.effect=='RR':
            if not self.algorithm in "MH,IV":
                Err="error.(algorithm for RR should be:'MH' or 'IV')"
                raise Exception(Err)
        elif self.effect=='RD':
            if not self.algorithm in "MH,IV":
                Err="error.(algorithm for RD should be:'MH' or 'IV')"
                raise Exception(Err)
        elif self.effect=='MD':
            if not self.algorithm in "IV":
                Err="error.(algorithm for MD should be:'IV')"
                raise Exception(Err)
        elif self.effect=='SMD':
            if not self.algorithm in "IV,IV-Heg,IV-Cnd,IV-Gls":
                Err="error.(algorithm for SMD should be:'IV','IV-Heg','IV-Cnd' or 'IV-Gls')"
                raise Exception(Err)

    #input:[[study],...], and the following attrs should be presetted before perform Meta()
    #models('Fixed','Random'),effect('OR','RR','RD','MD','SMD','SMD') and algorithm('MH','Peto','IV','IV-Heg','IV-Cnd','IV-Gls')
    #output:[[Total...],[study1...],...]
    def meta0 (self, stds):  
        self.chkmethods()
        if self.datatype=='CATE':
            self.studies=CATE_en2abcd(stds)  #so, Data.studies may diff to Meta.studies
        elif self.datatype=='CONT':
            self.studies=stds

        if self.models=='Fixed':  #Default
            if self.algorithm=='MH' and self.effect=='OR' :
                result=MH_total_OR (self.studies)
            elif self.algorithm=='MH' and self.effect=='RR' :
                result=MH_total_RR (self.studies)
            elif self.algorithm=='MH' and self.effect=='RD' :
                result=MH_total_RD (self.studies)
            elif self.algorithm=='Peto' and self.effect=='OR' :
                result=Peto_total_OR (self.studies)
            elif self.algorithm=='IV' and self.effect=='OR' :
                result=IV_total_OR (self.studies)
            elif self.algorithm=='IV' and self.effect=='RR' :
                result=IV_total_RR (self.studies)
            elif self.algorithm=='IV' and self.effect=='RD' :
                result=IV_total_RD (self.studies)
            elif self.algorithm=='IV' and self.effect=='RD' :
                result=IV_total_RD (self.studies)
            elif self.algorithm=='IV' and self.effect=='MD' :
                result=IV_total_MD (self.studies)
            elif self.algorithm=='IV' and self.effect=='SMD' : #SMD default for 'Heg'
                result=IV_total_SMD (self.studies)
            elif self.algorithm=='IV-Heg' and self.effect=='SMD' :
                result=IV_total_SMD (self.studies,'Fixed','Heg')
            elif self.algorithm=='IV-Cnd' and self.effect=='SMD' :
                result=IV_total_SMD (self.studies,'Fixed','Cnd')
            elif self.algorithm=='IV-Gls' and self.effect=='SMD' :
                result=IV_total_SMD (self.studies,'Fixed','Gls')
            else :
                Err="error.(fail to get algorithm and effect settings in Fixed model)"
                raise Exception(Err)

        elif self.models=='Random':
            if self.algorithm=='MH' and self.effect=='OR' :
                result=MH_total_OR (self.studies,'Random')
            elif self.algorithm=='MH' and self.effect=='RR' :
                result=MH_total_RR (self.studies,'Random')
            elif self.algorithm=='MH' and self.effect=='RD' :
                result=MH_total_RD (self.studies,'Random')
            elif self.algorithm=='IV' and self.effect=='OR' :
                result=IV_total_OR (self.studies,'Random')
            elif self.algorithm=='IV' and self.effect=='RR' :
                result=IV_total_RR (self.studies,'Random')
            elif self.algorithm=='IV' and self.effect=='RD' :
                result=IV_total_RD (self.studies,'Random')
            elif self.algorithm=='IV' and self.effect=='MD' :
                result=IV_total_MD (self.studies,'Random')
            elif self.algorithm=='IV' and self.effect=='SMD' :  #SMD default for 'Heg'
                result=IV_total_SMD (self.studies,'Random')
            elif self.algorithm=='IV-Heg' and self.effect=='SMD' :
                result=IV_total_SMD (self.studies,'Random','Heg')
            elif self.algorithm=='IV-Cnd' and self.effect=='SMD' :
                result=IV_total_SMD (self.studies,'Random','Cnd')
            elif self.algorithm=='IV-Gls' and self.effect=='SMD' :
                result=IV_total_SMD (self.studies,'Random','Gls')
            else :
                Err="error.(fail to get algorithm and effect settings in Random model)"
                raise Exception(Err)
        else:
            Err="error.(fail to get algorithm and effect settings)"
            raise Exception(Err)

        return result

    #meta and subgroup analysis
    #input:studies[[study],...]
    #      subgrp[num,[subgroupname,index of studies,... ]...]
    #output:results[[Total...],[study1...],[subgroup1,...],[studyn,...]...[subgroupk,...]]
    def meta (self, studies, nosubgrp=False):
        if nosubgrp==True:         #force to no subgroup
            results=self.meta0 (studies)
            return results

        subgrp=self.subgroup
        if subgrp[0]<1 or subgrp[0]!=len(subgrp)-1:  #no subgroup data
            results=self.meta0 (studies)
            return results

        stds=[];rult=[]
        results=[self.meta0 (studies)[0]]  #get total

        for i in range(1, len(subgrp)):
            stds.extend(studies[i] for i in subgrp[i][1:])
            rult=self.meta0 (stds)
            rult[0][0]="<sub>"+subgrp[i][0]  #I'm a subgroup
            rult.append(rult[0])
            del rult[0]
            results.extend(rult)
            stds=[]
        return results
    
#class of plot: forest，funnel, sensitivity    
class Fig:    
    def __init__(self,size=[6,6],dpi=80): #set figure
        self.size=size #inchs
        self.dpi=dpi   #default:80pts
        self.title="Meta-analysis Results (Pymeta.com)"  
        self.nototal=False
        
        #plt.rcParams['font.sans-serif']=['SimHei'] #Chinese compatible 
                                                    #Chinese font do not support unicode (u'\u00b2')
                                                    #in this case, fontdict={'family' : 'serif'} may be helpful
        plt.rcParams['axes.unicode_minus']=False    #show minus "-"
        logger.info("Load Fig().")
        
    #Method of Drawing ForestPlot
    #input list structure:es_w_ci
    #results[0]:total; [0]'Total', [1]OR, [2]weight, [3]LCI, [4]UCI, [5]N, [6]SE{Ln(OR)}, [7]Q, [8]p, [9]I2 ...
    #results[i]:studies; [0]study name, [1]effect size, [2]weight, [3]LCI, [4]UCI, [5]n,[6]quality ...
    def forest (self,rults):
        myfig=Fig_Forest (self.size,self.dpi,rults, self.title,self.nototal)
        return myfig

    #Method of Drawing FunnelPlot
    #input list structure:[[Effect Size,SE(ln(effs))],...]
    #effs[0][0:2]:total
    #effs[idx][0:2]:studies
    def funnel (self,effs):
        myfig = Fig_Funnel (self.size,self.dpi,effs)
        return myfig

#########################
#general functions
#########################
def midstr(content,startStr='',endStr='',inc=False): 
    if content=="":
        return ""
    if startStr+endStr=='':
        return content
    if startStr=='':
        startIndex = 0
    else:
        startIndex = content.find(startStr)

    if startIndex>=0:
        startIndex += len(startStr)
        if endStr=='':
            endIndex = len(content)
        else:
            endIndex = content.find(endStr,startIndex)
        if endIndex>=0:
            if inc:
                startIndex=startIndex-len(startStr) if startIndex-len(startStr)>0 else 0
                endIndex=endIndex+len(endStr) if endIndex+len(endStr)<len(content) else len(content)
            return content[startIndex:endIndex]
    return ""

#check&format datatype
def chkdtype(datatype):
    dt=datatype.upper().replace("DATA","").replace("-"," ").replace("_"," ").strip()
    dt=dt.split(" ")[0]
    a_cate=["CAT","BIN","CATE","DICH","BINARY","CATEGORICAL","COUNT","DICHOTOMOUS"]
    a_cont=["CON","CONT","CONTINUOUS"]

    if dt in a_cate:
        dt="CATE"
    elif dt in a_cont:
        dt="CONT"
    else:
        Err="error.(no datatype setting:'CATE', 'CONT')"
        raise Exception(Err)
    return dt

################################
'''
Statistics classes and functions
'''
class lmtbl_chisquare:
        seeds=[0.99999,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,
        0.09,0.08,0.07,0.06,0.05,0.04,0.03,0.02,0.01,0.008,
        0.006,0.004,0.002,0.001,0.0001,0.00001]
        data=[
                [0,0.015790774,0.064184755,0.148471862,0.2749959,0.454936425,
                0.708326304,1.07419517,1.642375062,2.705543971,2.874373729,
                3.064902115,3.283020525,3.537384875,3.841459149,4.217884752,
                4.70929242,5.411894595,6.634896712,7.033474337,7.550302698,
                8.283815075,9.549535733,10.82756622,15.13670526,19.51142097],
                [0.00002,0.210721031,0.446287103,0.713349888,1.021651261,
                1.386294376,1.832581484,2.407945609,3.218875825,4.605170186,
                4.815891217,5.051457289,5.318520074,5.626821434,5.991464547,
                6.43775165,7.013115795,7.824046011,9.210340372,9.656627475,
                10.23199162,11.04292184,12.4292162,13.81551056,18.42068074,23.02585093],
                [0.001122583,0.584374375,1.005174016,1.423652248,1.86916841,
                2.365973893,2.946166194,3.664870358,4.641627502,6.251388457,
                6.491457584,6.758692523,7.060314098,7.406879995,7.814727764,
                8.311170822,8.947287448,9.837409286,11.34486668,11.82697381,
                12.44655101,13.31640864,14.79551705,16.26623615,21.10751346,25.90174974],
                [0.008957633,1.063623219,1.648776626,2.194698434,2.7528427,
                3.356694001,4.044626491,4.878432967,5.988616694,7.77944034,
                8.043435288,8.336531703,8.666428235,9.044368369,9.487729037,
                10.02551929,10.71189829,11.6678434,13.27670414,13.78895432,
                14.4458427,15.36561125,16.9237582,18.46682695,23.51274245,28.47325542],
                [0.032484399,1.610307989,2.342534315,2.999908149,3.655499646,
                4.351460222,5.131867373,6.064430245,7.289276183,9.236356938,
                9.52107416,9.836591308,10.19102793,10.59623215,11.07049775,
                11.64433189,12.37461851,13.38822261,15.08627247,15.62512207,
                16.31499161,17.27897694,18.90737739,20.51500566,25.74483196,30.85618994],
                [0.07907433,2.20413068,3.070088414,3.827551604,4.570153831,
                5.348120843,6.210757195,7.231135332,8.558059721,10.64464068,
                10.94790172,11.28349557,11.65992262,12.08957827,12.59158724,
                13.19781465,13.96761693,15.03320775,16.81189383,17.37484553,
                18.09463447,19.09879292,20.79116772,22.45774449,27.85634124,33.10705682],
                [0.152859753,2.833106932,3.822321961,4.67133054,5.493234987,
                6.345811373,7.283207505,8.38343064,9.803249854,12.01703656,
                12.33723558,12.69117599,13.0877088,13.53973354,14.06714043,
                14.70304666,15.5090897,16.62242187,18.47530691,19.06047255,
                19.80786049,20.84911787,22.60067086,24.32188634,29.87750391,35.25853642],
                [0.255357772,3.489539134,4.593573645,5.527422145,6.422645649,
                7.344121629,8.350525468,9.524458194,11.03009143,13.36156614,
                13.69745554,14.06839705,14.4835684,14.95633906,15.50731306,
                16.17077561,17.01049321,18.16823077,20.09023503,20.69611949,
                21.46926575,22.54517756,24.35208135,26.12448156,31.827628,37.33159364],
                [0.386458863,4.168159042,5.380053335,6.393306001,7.35703456,
                8.342832783,9.413640132,10.65637209,12.24214549,14.68365662,
                15.0342307,15.42108788,15.85371313,16.34591786,16.91897762,
                17.60827685,18.47958642,19.6790161,21.66599433,22.29137421,
                23.08877044,24.19732982,26.05643335,27.87716488,33.71994844,39.34065373],
                [0.545169541,4.865182068,6.17907932,7.267218288,8.295471948,
                9.34181805,10.47323623,11.78072263,13.44195758,15.98717917,
                16.35160466,16.75347765,17.20257397,17.71312357,18.30703805,
                19.02074335,19.92191001,21.16076754,23.20925116,23.8531005,
                24.67348033,25.81299977,27.72164723,29.58829845,35.56401394,41.29615797],
                [0.730073364,5.57778484,6.988673545,8.147867844,9.237285532,
                10.34099825,11.52983382,12.89866814,14.63142049,17.27500852,
                17.65258016,18.06870697,18.53344169,19.06141276,19.67513757,
                20.41203411,21.34158304,22.6179408,24.72497031,25.38641214,
                26.22869205,27.39772587,29.35363761,31.26413362,37.36698644,43.20595972],
                [0.939594769,6.303796082,7.807327768,9.034276765,10.18197166,
                11.34032282,12.58383797,14.01110017,15.81198622,18.54934779,
                18.93945857,19.36918308,19.84883885,20.39343548,21.02606982,
                21.78510904,22.74176613,24.05395669,26.21696731,26.89524046,
                27.75847919,28.95577019,30.95696053,32.90949041,39.13440388,45.07614652],
                [1.17214265,7.041504637,8.633860877,9.925682506,11.12914009,
                12.33975614,13.63557102,15.11872166,16.98479702,19.81192931,
                20.21404951,20.65679894,21.1507493,21.71127608,22.3620325,
                23.14229721,24.12494701,25.47150915,27.68824961,28.38268506,
                29.26605478,30.49049616,32.53521451,34.52817898,40.87065501,46.91155313],
                [1.426184456,7.789533743,9.467328089,10.82147793,12.07848283,
                13.33927471,14.68529426,16.22209861,18.15077056,21.06414421,
                21.47780626,21.93307486,22.44076648,23.01660895,23.68479131,
                24.48547022,25.49312548,26.87276464,29.14123774,29.85124128,
                30.75400624,32.00461359,34.09130099,36.12327368,42.57928895,48.7160969],
                [1.700281207,8.546756297,10.30695923,11.72116908,13.02974978,
                14.33885982,15.73322295,17.32169449,19.31065711,22.30712958,
                22.73191673,23.19925365,23.72019308,24.31080308,24.99579013,
                25.81615891,26.84793759,28.25949634,30.57791417,31.30295333,
                32.22445389,33.50034396,35.62760013,37.69729823,44.26322495,50.49300558],
                [1.993100942,9.312236471,11.15211657,12.62434898,13.98273671,
                15.33849951,16.77953671,18.41789439,20.46507929,23.54182892,
                23.97736669,24.45636665,24.99011023,25.59499477,26.29622761,
                27.13563426,28.19074206,29.63317731,31.99992691,32.73952089,
                33.67916114,34.97953551,37.14609345,39.25235481,45.92489905,52.24497689],
                [2.303420576,10.08518638,12.00226593,13.53067656,14.93727199,
                16.33818271,17.82438727,19.51102236,21.61456054,24.76903535,
                25.21498476,25.70528088,26.25142689,26.87013943,27.58711164,
                28.44496525,29.52268192,30.99504719,33.4086636,34.16237538,
                35.11961341,36.44374562,38.64845152,40.79021671,47.56636956,53.97429344],
                [2.630121722,10.86493621,12.85695319,14.43986255,15.89321209,
                17.337903,18.86790412,20.60135412,22.75954582,25.9894231,
                26.44547544,26.9467335,27.5049159,28.13704951,28.86929943,
                29.74506112,30.84472954,32.34616093,34.80530572,35.57273575,
                36.54707593,37.89430122,40.13609844,42.31239633,49.18939448,55.68290738],
                [2.972183685,11.65091021,13.71578988,15.35166063,16.85043363,
                18.33765323,19.91019886,21.68912657,23.90041722,27.20357106,
                27.6694436,28.18135749,28.75124092,29.3964227,30.14352721,
                31.03670291,32.15772038,33.68742507,36.19086911,36.97165019,
                37.96263691,39.33234389,41.61025992,43.82019596,50.79548968,57.37250401],
                [3.328675327,12.44260928,14.57843952,16.26585666,17.80882981,
                19.33742983,20.95136838,22.77454507,25.03750563,28.41198058,
                28.88741314,29.40970112,29.99097688,30.64886346,31.41043286,
                32.3205674,33.46237849,35.01962554,37.56623475,38.36002791,
                39.36724042,40.75886424,43.07200006,45.31474662,52.38597329,59.04455039],
                [3.698746821,13.2395981,15.44460854,17.1822655,18.76830962,
                20.33722858,21.9914975,23.8577889,26.17109991,29.61508943,
                30.09984149,30.63224298,31.22462565,31.89489954,32.67057337,
                33.59724563,34.75933632,36.34344894,38.93217269,39.73866365,
                40.76171179,42.17472852,44.52224969,46.79703804,53.96200014,60.70033312],
                [4.081621648,14.04149341,16.31404003,18.10072389,19.7287913,
                21.33704534,23.03066093,24.93901574,27.30145403,30.81328234,
                31.3071307,31.84940398,32.45262859,33.13499492,33.92443852,
                34.86725738,36.04915009,37.65949929,40.28936044,41.10825694,
                42.14677757,43.58069941,45.96182867,48.26794229,55.5245888,62.3409881],
                [4.476589238,14.84795589,17.18650624,19.02108741,20.69020469,
                22.3368793,24.06892481,26.01836515,28.42879253,32.00689967,
                32.50963648,33.06155637,33.67537609,34.36956026,35.17246163,
                36.13106262,37.33231213,38.9683113,41.63839812,42.46942719,
                43.52308136,44.97745264,47.39146361,49.72823246,57.07464314,63.9675242],
                [4.882998242,15.6586842,18.0618045,19.94322915,21.65248637,
                23.33672677,25.10634822,27.09596132,29.55331525,33.19624426,
                33.70767532,34.2690314,34.89321559,35.59896088,36.4150285,
                37.38907067,38.60926067,40.27036103,42.97982015,43.82272622,
                44.89119634,46.36559023,48.81180207,51.17859777,58.61296975,65.58084238],
                [5.300250456,16.47340823,18.93975474,20.86703465,22.61557929,
                24.33658743,26.14298396,28.17191525,30.67520091,34.38158698,
                34.90153029,35.47212537,36.10645783,36.82352361,37.65248413,
                38.64164792,39.88038783,41.56607452,44.31410491,45.16864801,
                46.25163585,47.74565138,50.22342413,52.61965576,60.14029191,67.18175124],
                [5.727795481,17.29188535,19.82019439,21.79240089,23.57943433,
                25.33645925,27.17887957,29.24632692,31.7946101,35.56317121,
                36.09145577,36.67110455,37.31538211,38.0435422,38.88513865,
                39.88912356,41.14604605,42.85583483,45.64168268,46.5076369,
                47.60486163,49.11812107,51.62685184,54.05196237,61.65726128,68.77097981],
                [6.165125825,18.11389612,20.70297662,22.71923661,24.54400558,
                26.3363399,28.21407802,30.31928641,32.91168775,36.74121675,
                37.27768139,37.86620931,38.52024053,39.25928183,40.11327205,
                41.13179486,42.40655379,44.13998791,46.96294214,47.84009429,
                48.95129088,50.48343765,53.02255705,55.47602018,63.16446742,70.34918809],
                [6.61177268,18.93924261,21.58796957,23.64745785,25.50925112,
                27.3362301,29.24861825,31.39087543,34.0265652,37.91592255,
                38.46041536,39.05765747,39.72126161,40.47098284,41.33713813,
                42.36993129,43.66219955,45.41884745,48.27823579,49.16638427,
                50.29130203,51.8419988,54.410968,56.89228537,64.66244582,71.91697595],
                [7.067302126,19.76774391,22.47505242,24.57698803,26.47513514,
                28.3361282,30.28253597,32.4611683,35.1393618,39.08746978,
                39.63984679,40.24564731,40.9186532,41.6788639,42.55696777,
                43.60377803,44.91324593,46.6926988,49.5878845,50.48683826,
                51.62523959,53.19416659,55.79247454,58.30117346,66.15168462,73.47489052],
                [7.531311771,20.59923476,23.36411522,25.50775906,27.4416231,
                29.33603221,31.31586326,33.53023285,36.25018677,40.25602376,
                40.81614856,41.43035953,42.11260506,42.88312462,43.77297178,
                44.83355892,46.15993268,47.96180281,50.89218135,51.80175901,
                52.95341831,54.54027183,57.16743307,59.70306426,67.63263025,75.02343247],
                [8.003427729,21.43356472,24.25505672,26.43970716,28.40868306,
                30.33594347,32.34862996,34.59813124,37.35913986,41.42173586,
                41.98947889,42.61195987,43.30329113,44.08394782,44.98534322,
                46.059479,47.40247944,49.22639823,52.19139488,53.11142395,
                54.27612662,55.88061772,58.53617025,61.09830603,69.10569228,76.56306137],
                [8.4833021,22.27059479,25.14778555,27.37277355,29.37628808,
                31.33586054,33.38086336,35.66492071,38.46631278,42.58474512,
                43.15998312,43.79060049,44.49087077,45.28150155,46.19425944,
                47.28172663,48.64108798,50.48670449,53.4857719,54.41608805,
                55.59362964,57.21548299,59.89898645,62.487219,70.57124757,78.09420036],
                [8.970610441,23.11019718,26.04221622,28.30690588,30.34441134,
                32.33578165,34.41258863,36.73065353,39.57178996,43.74517961,
                44.32779512,44.96642152,45.67549112,46.47594027,47.39988381,
                48.50047532,49.8759442,51.74292377,54.77553984,55.71598628,
                56.90617171,58.54512455,61.25615857,63.87009853,72.02964378,79.61723994],
                [9.465049772,23.95225346,26.93826965,29.24205434,31.31302874,
                33.33570845,35.443829,37.79537838,40.67564938,44.90315759,
                45.49303851,46.1395524,46.85728797,47.66740695,48.60236738,
                49.71588534,51.10721981,52.99524285,56.06090884,57.01133579,
                58.21397862,59.86977982,62.60794255,65.24721747,73.4812025,81.13254149],
                [9.966336672,24.79665505,27.83587413,30.17817328,32.28211671,
                34.33563964,36.47460575,38.85914021,41.77796322,46.05878853,
                46.65582772,47.31011294,48.03638703,48.85603384,49.80184958,
                50.92810508,52.33507381,54.24383469,57.34207343,58.30233767,
                59.51725954,61.18966877,63.95457553,66.61882885,74.92622188,82.64044019],
                [10.47420561,25.64330025,28.73496185,31.11521954,33.2516567,
                35.3355749,37.50493927,39.92198072,42.87879847,47.21217401,
                47.81626896,48.47821432,49.21290495,50.04194369,50.99846018,
                52.13727227,53.55965378,55.48885983,58.61921449,59.58917865,
                60.81620865,62.50499562,65.29627774,67.98516763,76.36497895,84.14124764],
                [10.98840738,26.49209476,29.63547018,32.05315524,34.22162857,
                36.33551222,38.53484797,40.98393868,43.97821737,48.3634085,
                48.97446096,49.64395994,50.38695021,51.22525064,52.19231975,
                53.34351515,54.78109696,56.73046761,59.89250003,60.8720323,
                62.11100663,63.81595041,66.63325413,69.34645251,77.79773169,85.63525411],
                [11.50870795,27.34295027,30.53734011,32.99194311,35.1920133,
                37.33545422,39.56434917,42.04505023,45.07627797,49.51257981,
                50.13049574,50.80744614,51.5586239,52.40606104,53.38354065,
                54.54695278,55.99953128,57.96879722,61.16208675,62.15106079,
                63.40182191,65.12271033,67.96569581,70.70288742,79.2247208,87.12273053],
                [12.03488698,28.19578549,31.44051828,33.93154899,36.16279619,
                38.33539938,40.59345905,43.10534911,46.1730347,50.65977047,
                51.2844592,51.96876287,52.72802039,53.58447423,54.5722278,
                55.74769677,57.21507619,59.20397864,62.428121,63.42641527,
                64.68881186,66.4254409,69.29378135,72.05466296,80.64617127,88.60393026],
                [12.56673691,29.05052334,32.34495328,34.87193994,37.13396078,
                39.33534592,41.62219275,44.16486689,47.26853774,51.80505719,
                52.43643149,53.12799426,53.89522798,54.76058309,55.75847932,
                56.94585138,58.4278436,60.43613349,63.69073973,64.69823728,
                65.97212358,67.72429702,70.61767787,73.40195753,82.06229381,90.07909063],
                [13.1040619,29.90709191,33.2505974,35.81308776,38.10549251,
                40.33529635,42.65056447,45.22363278,48.36283478,52.94851197,
                53.58648828,54.28521916,55.06032939,55.93447468,56.9423872,
                58.14151432,59.63793796,61.66537576,64.95007131,65.96665949,
                67.25189533,69.01942393,71.93754207,74.74493841,83.47328607,91.54843436],
                [13.64667698,30.76542371,34.15740519,36.75496361,39.07737635,
                41.33524929,43.67858751,46.28167511,49.45597036,54.09020241,
                54.73470006,55.44051131,56.22340226,57.1062307,58.12403775,
                59.33477731,60.8454577,62.89181243,66.20623626,67.23180648,
                68.52825661,70.31095798,73.2535211,76.08376273,84.87933374,93.01217076],
                [14.19440719,31.62545472,35.06533579,37.69754113,40.04960178,
                42.3352046,44.70627478,47.33901971,50.54798635,55.23019204,
                55.88113324,56.5939407,57.38451957,58.27592795,59.3035121,
                60.52572656,62.05049525,64.11554411,67.45934789,68.49379537,
                69.80132928,71.59902742,74.56575338,77.41857826,86.28061157,94.47049685],
                [14.74708683,32.48712621,35.97434906,38.64079443,41.02215568,
                43.33516031,45.73363743,48.39569113,51.63892217,56.36854066,
                57.02585036,57.7455726,58.5437498,59.44363875,60.48088667,
                61.71444326,63.25313762,65.33666553,68.70951293,69.75273648,
                71.07122813,72.88375303,75.87436927,78.74952425,87.67728423,95.92359836],
                [15.30455901,33.35038143,36.88440769,39.58470202,41.99502667,
                44.33511948,46.76068686,49.45171252,52.72881496,57.50530467,
                58.16891036,58.89546875,59.70115815,60.60943126,61.65623348,
                62.90100395,64.45346688,66.55526601,69.95683202,71.00873385,
                72.33806141,74.1652486,77.17949176,80.07673204,89.06950709,97.37165058],
                [15.8666747,34.2151672,37.79547647,40.52924097,42.96820255,
                45.33508057,47.78743358,50.50710572,53.81769981,58.64053728,
                59.3103689,60.04368742,60.85680559,61.77336985,62.82962054,
                64.08548088,65.65156049,67.77143007,71.2014002,72.26188571,
                73.60193134,75.44362199,78.481237,81.40032569,90.45742699,98.81481923],
                [16.43329239,35.08143274,38.70752123,41.47439045,43.94167575,
                46.33504186,48.81388751,51.5618914,54.90560987,59.77428882,
                60.45027858,61.19028371,62.01074991,62.93551511,64.00111212,
                65.26794238,66.84749168,68.9852371,72.44330732,73.51228497,
                74.8629346,76.71897481,79.77971485,82.72042255,91.84118284,100.2532611],
                [17.00427753,35.94913143,39.62051203,42.42012961,44.91543551,
                47.33500609,49.84005805,52.61608909,55.99257651,60.90660689,
                61.58868918,62.33530973,63.16304584,64.09592508,65.17076907,
                66.4484531,68.0413298,70.19676276,73.68263846,74.76001955,
                76.1211627,77.99140339,81.07502931,84.03713375,93.22090622,101.6871247],
                [17.57950209,36.81821779,40.53441833,43.36644173,45.88947274,
                48.3349719,50.86595411,53.66971732,57.07862945,62.03753663,
                62.72564788,63.4788149,64.31374522,65.25465411,66.33864905,
                67.62707431,69.23314052,71.40607857,74.91947424,76.00517284,
                77.3767024,79.26099906,82.36727901,85.35056465,94.59672192,103.1165511],
                [18.15884409,37.68864904,41.44921183,44.31330821,46.86377731,
                49.33493921,51.89158414,54.72279365,58.16379654,63.16712082,
                63.86119943,64.62084608,65.46289724,66.41175397,67.50480652,
                68.80386416,70.42298619,72.61325241,76.15389117,77.24782391,
                78.62963601,80.52784846,83.65655752,86.66081523,95.96874841,104.5416741],
                [18.74218722,38.56038459,42.36486481,45.26071223,47.83834375,
                50.33490621,52.91695616,55.77533522,59.2481052,64.29540012,
                64.99538633,65.76144778,66.61054859,67.56727382,68.66929388,
                69.97887784,71.61092601,73.81834874,77.38596193,78.48804792,
                79.88004174,81.79203392,84.94295377,87.96798052,97.33709828,105.9626208],
                [19.32942044,39.43338524,43.28135324,46.20863668,48.81316305,
                51.33487592,53.94207733,56.82735709,60.33158059,65.42241317,
                66.12824898,66.90066231,67.75674369,68.7212604,69.83216031,
                71.15216786,72.79701625,75.02142887,78.61575562,79.7259163,
                81.12799395,83.05363374,86.22655232,89.27215088,98.70187866,107.3795123],
                [19.92043793,40.30761533,44.19865214,47.15706833,49.78822784,
                52.33484687,54.96695574,57.87887472,61.41424693,66.54819705,
                67.25982583,68.03852993,68.90152479,69.87375818,70.99345279,
                72.32378416,73.9813105,76.22255111,79.84333801,80.96149706,
                82.37356343,84.31272243,87.50743366,90.57341236,100.0631916,108.792464],
                [20.51513831,41.18303937,45.11673836,48.10599181,50.76352953,
                53.33481902,55.99159793,58.92990249,62.49612729,67.6727862,
                68.39015349,69.175089,70.04493215,71.02480953,72.15321612,
                73.49377434,75.16385977,77.42177107,81.06877179,82.19485502,
                83.61681762,85.56937102,88.78567438,91.87184694,101.4211343,110.2015857],
                [21.11342469,42.05962401,46.03558983,49.0553933,51.73906417,
                54.33479232,57.01601051,59.98045409,63.57724367,68.79621429,
                69.5192672,70.31037608,71.18700415,72.17445485,73.31149298,
                74.6621835,76.34471273,78.61914172,82.2921167,83.42605195,
                84.85782084,86.82364721,90.06134791,93.16753284,102.7757994,111.6069825],
                [21.71520431,42.93733742,46.95518456,50.00525963,52.71482412,
                55.33476454,58.04019975,61.03054262,64.65761705,69.91851319,
                70.64719969,71.44442606,72.3277774,73.32273268,74.4683241,
                75.82905545,77.52391582,79.81471365,83.51342979,84.65514683,
                86.09663449,88.07561567,91.33452397,94.46054471,104.1272756,113.0087544],
                [22.32038832,43.81614853,47.87550416,50.95557695,53.69080326,
                56.33473956,59.06417168,62.08018056,65.73726747,71.03971332,
                71.77398307,72.57727228,73.46728689,74.46967985,75.6237484,
                76.9944313,78.70151339,81.00853518,84.73276555,85.88219602,
                87.33331722,89.32533814,92.6052693,95.75095391,105.4756476,114.4069973],
                [22.92889154,44.69602945,48.79652863,51.90633555,54.66699425,
                57.33471556,60.08793204,63.12937986,66.8162141,72.15984378,
                72.89964774,73.70894691,74.60556605,75.61533156,76.77780308,
                78.15835061,79.87754785,82.2006525,85.95017607,87.10725339,
                88.56792511,90.57287369,93.8736477,97.03882865,106.8209962,115.8018028],
                [23.54063226,45.57695184,49.71823959,52.85752288,55.64339423,
                58.33469248,61.11148633,64.17815192,67.89447526,73.2789324,
                74.02422281,74.83947983,75.74264687,76.75972147,77.93052372,
                79.32085118,81.05205977,83.39110979,87.16571142,88.33037047,
                89.80051182,91.81827882,95.13972024,98.32423422,108.1633991,117.1932586],
                [24.15553194,46.45888908,50.64061945,53.80912793,56.61999662,
                59.33466815,62.13483981,65.2265077,68.97206851,74.39700583,
                75.14773614,75.96890051,76.87855999,77.90288182,79.08194439,
                80.48196925,82.22508799,84.57994937,88.37941893,89.55159664,
                91.03112874,93.06160765,96.40354541,99.60723316,109.5029305,118.5814487],
                [24.77351543,47.34181565,51.56365033,54.76114016,57.59679629,
                60.33464644,63.15799752,66.27445767,70.04901064,75.51408958,
                76.27021442,77.09723689,78.01333473,79.0448435,80.23209774,
                81.64173952,83.39666972,85.7672118,89.59134452,90.77097919,
                92.25982511,94.30291204,97.66517926,100.8878854,110.8396615,119.9664538],
                [25.39451023,48.22570711,52.48731794,55.71354805,58.57378835,
                61.33462553,64.1809643,67.32201187,71.12531778,76.63020814,
                77.39168324,78.22451578,79.14699958,80.18563612,81.38101507,
                82.80019525,84.56684039,86.95293594,90.80153206,91.98856348,
                93.48664815,95.54224172,98.92467554,102.1662484,112.1736605,121.3483511],
                [26.0184467,49.11053923,53.41160592,56.6663445,59.55096638,
                62.33460539,65.20374543,68.36917939,72.20100539,77.74538499,
                78.51216716,79.35076291,80.27958089,81.32528811,82.52872641,
                83.95736837,85.73563476,88.13715911,92.01002365,93.20439304,
                94.71164318,96.77964443,100.1820859,103.4423774,113.5049929,122.7272149],
                [26.64525782,49.99629088,54.33649938,57.619519,60.52832899,
                63.3345839,66.22634421,69.41597049,73.27608878,78.85964267,
                79.63168975,80.47600295,81.41110487,82.46382674,83.6752606,
                85.11328954,86.90308547,89.31991713,93.2168597,94.4185099,
                95.93485373,98.01516601,101.4374597,104.7163254,114.8337217,124.1031165],
                [27.27487907,50.88293976,55.26198404,58.57306262,61.5058704,
                64.33456486,67.24876545,70.46239359,74.35058134,79.97300285,
                80.75027367,81.60025965,82.54159644,83.60127823,84.82064534,
                86.26798822,88.0692241,90.50124442,94.42207905,95.6309538,
                97.15632161,99.2488505,102.6908448,105.9881432,116.1599075,125.4761244],
                [27.9072483,51.77046507,56.18804509,59.52696677,62.48358645,
                65.33454647,68.27101327,71.50845724,75.4244972,81.08548634,
                81.86794068,82.72355579,83.67107952,84.73766777,85.96490727,
                87.42149271,89.23408104,91.68117408,95.62571904,96.84176355,
                98.37608705,100.4807403,103.9422869,107.2578799,117.4836084,126.8463046],
                [28.54230559,52.6588468,57.11467126,60.48122169,63.46147313,
                66.33452873,69.29309164,72.55416966,76.49784953,82.19711318,
                82.98471174,83.84591332,84.79957715,85.87301998,87.108072,
                88.57383026,90.39768549,92.85973794,96.82781561,98.05097633,
                99.59418871,101.710876,105.19183,108.5255824,118.8048806,128.2137207],
                [29.17999316,53.54806568,58.04184893,61.43582218,64.43952472,
                67.33451161,70.31500435,73.59953875,77.57065104,83.30790264,
                84.100607,84.96735336,85.92711149,87.00735749,88.25016421,
                89.72502708,91.56006559,94.03696665,98.02840334,99.25862789,
                100.8106638,102.9392971,106.4395166,109.7912965,120.1237779,129.5784337],
                [29.82025513,54.43810229,58.96956585,62.39075931,65.41774094,
                68.33449278,71.33675508,74.64457214,78.64291394,84.4178733,
                85.21564588,86.08789624,87.05370388,88.1407031,89.39120764,
                90.87510843,92.72124842,95.21288972,99.22751553,100.4647527,
                102.0255483,104.1660412,107.6853875,111.0550655,121.4403522,130.9405025],
                [30.46303792,55.32894035,59.89781021,63.34602571,66.39611662,
                69.33447649,72.35834738,75.68927715,79.71465002,85.52704303,
                86.32984708,87.20756155,88.17937489,89.27307847,90.53122518,
                92.02409864,93.88126009,96.38753561,100.4251843,101.6693839,
                103.2388768,105.3911447,108.9294821,112.3169318,122.7546538,132.2999839],
                [31.10828945,56.22056189,60.82657064,64.30161426,67.37464831,
                70.33446073,73.37978464,76.73366088,80.78587063,86.6354291,
                87.44322865,88.32636818,89.30414436,90.40450447,91.6702389,
                93.17202117,95.04012576,97.56093173,101.6214406,102.8725536,
                104.4506827,106.6146428,110.1718384,113.5769359,124.0667309,133.6569329],
                [31.7559595,57.11295026,61.75583495,65.25751657,68.3533327,
                71.3344455,74.40107017,77.77773015,81.85658672,87.74304768,
                88.55580797,89.44433437,90.42803141,91.53500124,92.80827009,
                94.31889866,96.19786972,98.73310454,102.8163143,104.0742927,
                105.6609978,107.8365694,111.412493,114.8351171,125.3766302,135.0114021],
                [32.4059995,58.00608943,62.68559481,66.21372886,69.3321646,
                72.33443079,75.42220714,78.82149154,82.92680886,88.8499157,
                89.66760184,90.5614777,91.55105451,92.66458822,93.94533966,
                95.46475294,97.35451539,99.90407929,104.0098342,105.2746309,
                106.8698534,109.0569572,112.6514813,116.0915131,126.6843966,136.3634427],
                [33.0583625,58.89996391,63.61583892,67.17024341,70.31114456,
                73.33441409,76.4431986,79.86495141,83.99654726,89.95604824,
                90.77862649,91.67781514,92.67323148,93.79328418,95.08146673,
                96.60960511,98.51008542,101.0738812,105.2020281,106.4735971,
                108.0772795,110.2758377,113.8888375,117.3461605,127.9900738,137.7131041],
                [33.71300307,59.7945578,64.54655743,68.12705407,71.29026789,
                74.3344,77.46404754,80.90811592,85.06581178,91.06146028,
                91.88889756,92.79336312,93.79457954,94.92110725,96.21667082,
                97.75347556,99.66460165,102.2425339,106.392923,107.671219,
                109.283305,111.4932416,115.1245945,118.5990947,129.2937036,139.0604339],
                [34.36987722,60.68985835,65.47774078,69.08415486,72.26953169,
                75.33438635,78.48475595,81.950991,86.13461195,92.16616631,
                92.99842974,93.90813748,94.91511534,96.04807495,97.35097045,
                98.89638396,100.8180852,103.4100602,107.5825449,108.8675237,
                110.4879579,112.7091982,116.3587843,119.8503498,130.5953269,140.4054782],
                [35.02894239,61.58585076,66.40937856,70.04153848,73.24893319,
                76.33437315,79.50532819,82.99358239,87.20295699,93.27018032,
                94.10723854,95.02215355,96.03485498,97.17420424,98.48438354,
                100.0383494,101.9705566,104.5764826,108.7709188,110.062537,
                111.6912653,113.9237362,117.5914378,121.0999588,131.8949828,141.7482816],
                [35.69015729,62.4825215,67.34146404,70.99920221,74.22846766,
                77.33436038,80.52576619,84.03589641,88.2708558,94.37351585,
                95.21533772,96.13542618,97.15381403,98.29951149,99.61692741,
                101.1793902,103.1220355,105.7418226,109.9580692,111.2562841,
                112.8932535,115.1368831,118.8225848,122.3479538,133.1927095,143.0888873],
                [36.35348179,63.37985748,68.27398754,71.95713931,75.20813622,
                78.33434542,81.54607256,85.07793701,89.33831701,95.476186,
                96.32274096,97.24796972,98.27200756,99.42401259,100.7486188,
                102.3195244,104.2725411,106.906101,111.1440195,112.4487895,
                114.0939478,116.3486657,120.0522541,123.5943655,134.4885437,144.4273372],
                [37.01887744,64.27784601,69.20694071,72.91534455,76.1879346,
                79.33433312,82.56624978,86.11971006,90.40534898,96.57820347,
                97.42946155,98.35979807,99.38945015,100.5477229,101.8794741,
                103.458769,105.4220919,108.069338,112.3287926,113.6400766,
                115.2933727,117.5591098,121.280474,124.839224,135.7825211,145.7636716],
                [37.6863065,65.17647386,70.14031544,73.87381289,77.16786035,
                80.3343212,83.58630031,87.16122062,91.4719598,97.67958055,
                98.53551233,99.47092422,100.506156,101.6706573,103.0095088,
                104.5971414,106.5707058,109.2315532,113.5124106,114.8301683,
                116.4915521,118.7682406,122.5072713,126.0825583,137.0746764,147.0979299],
                [38.35573248,66.0757309,71.07410389,74.83253941,78.14791108,
                81.33430965,84.60622648,88.20247356,92.53815665,98.78032915,
                99.64090576,100.5813621,101.6221387,102.7928302,104.1387383,
                105.7346568,107.7184003,110.3927655,114.6948948,116.0190868,
                117.6885089,119.9760825,123.7326725,127.3243966,138.365043,148.4301501],
                [39.02711998,66.97560475,72.00829716,75.79151766,79.1280825,
                82.33429846,85.62603059,89.24347363,93.60394836,99.88046081,
                100.7456539,101.691124,102.7374115,103.9142555,105.2671774,
                106.8713314,108.8651922,111.5529934,115.876266,117.2068536,
                118.8842657,121.1826587,124.956703,128.564766,139.6536534,149.7603691],
                [39.70043464,67.87608423,72.94289026,76.75074621,80.10837611,
                83.33428497,86.64571485,90.28422542,94.66934173,100.9799867,
                101.8497685,102.8002221,103.8519875,105.0349469,106.3948404,
                108.0071805,110.0110978,112.7122545,117.0565444,118.3934894,
                120.078844,122.3879922,126.1793875,129.8036932,140.9405391,151.0886226],
                [40.37564307,68.77715848,73.87787491,77.7102191,81.08878799,
                84.33427414,87.66528141,91.32473337,95.73434389,102.0789178,
                102.9532609,103.9086683,104.965879,106.1549175,107.5217411,
                109.1422187,111.1561331,113.8705663,118.2357494,119.5790145,
                121.2722648,123.5921052,127.40075,131.0412037,142.2257307,152.4149455],
                [41.05271284,69.67881696,74.81324416,78.66993199,82.06931603,
                85.33426363,88.68473234,92.36500182,96.79896178,103.1772645,
                104.0561421,105.0164742,106.0790981,107.2741801,108.6478931,
                110.2764606,112.3003134,115.0279455,119.4139,120.7634484,
                122.4645487,124.7950192,128.6208138,132.2773225,143.509258,153.7393714],
                [41.73161242,70.58104946,75.74899128,79.62988066,83.04995821,
                86.33425344,89.70406969,93.40503497,97.86320214,104.2750371,
                105.1584228,106.123651,107.1916562,108.392747,109.7733095,
                111.4099199,113.4436537,116.1844084,120.5910147,121.9468102,
                123.6557154,125.9967551,129.8396016,133.5120737,144.7911497,155.061933],
                [42.41231117,71.48384496,76.6851097,80.59006103,84.03071254,
                87.33424356,90.7232954,94.44483689,98.92707151,105.3722455,
                106.2601134,107.2302096,108.3035656,109.5106303,110.898003,
                112.5426102,114.5861684,117.3399708,121.7671112,123.1291184,
                124.8457842,127.1973332,131.0571351,134.745481,146.0714338,156.3826623],
                [43.09477908,72.38719583,77.62159308,81.5504673,85.01157488,
                88.33423399,91.7424114,95.48441154,99.99057624,106.4688994,
                107.3612241,108.3361605,109.414837,110.6278416,112.0219859,
                113.6745447,115.7278717,118.4946481,122.942207,124.310391,
                126.0347737,128.3967733,132.2734359,135.977567,147.3501376,157.70159],
                [43.77898755,73.29109173,78.55843375,82.51109904,85.99254753,
                89.33422171,92.76142054,96.52376276,101.0537225,107.5650082,
                108.4617645,109.441514,110.5254811,111.7443924,113.1452703,
                114.8057361,116.8687773,119.6484552,124.1163189,125.4906453,
                127.2227023,129.5950944,133.4885246,137.2083541,148.6272876,159.0187465],
                [44.46490818,74.19552356,79.49562848,83.47195088,86.97362669,
                90.33421238,93.78032273,97.56289431,102.1165163,108.660581,
                109.5617444,110.5462802,111.6355085,112.8602936,114.2678679,
                115.9361969,118.0088985,120.8014067,125.2894633,126.6698984,
                128.4095875,130.7923153,134.7024214,138.4378637,149.9029094,160.3341607],
                [45.15251352,75.10048244,80.43317011,84.43301916,87.95481062,
                91.33420332,94.79912062,98.60180982,103.1789635,109.7556267,
                110.6611728,111.6504687,112.7449293,113.9755561,115.3897899,
                117.0659391,119.1482483,121.9535167,126.4616562,127.8481669,
                129.5954466,131.988454,135.9151459,139.666117,151.1770281,161.6478611],
                [45.84177683,76.00595976,81.37105297,85.39430035,88.93609758,
                92.33419452,95.81781593,99.64051282,104.2410697,110.8501538,
                111.7600589,112.754089,113.8537535,115.0901902,116.5110475,
                118.1949744,120.2868392,123.1047989,127.6329131,129.0254666,
                130.7802962,133.1835281,137.1267172,140.8931342,152.4496681,162.9598755],
                [46.53267208,76.9119461,82.30927151,86.355791,89.91748592,
                93.33418599,96.83641031,100.6790059,105.3028405,111.9441708,
                112.8584115,113.8571503,114.9619907,116.2042062,117.6316514,
                119.3233144,121.424684,124.2552667,128.8032491,130.2018133,
                131.9641527,134.3775549,138.3371538,142.1189353,153.720853,164.2702307],
                [47.22517389,77.81843523,83.24782037,87.31748581,90.89897162,
                94.33417479,97.85490538,101.7172941,106.3642812,113.0376858,
                113.956239,114.9596615,116.0696504,117.3176133,118.751612,
                120.4509701,122.5617938,125.4049332,129.972679,131.3772222,
                133.1470318,135.5705509,139.5464738,143.3435397,154.9906059,165.578953],
                [47.91925756,78.72541827,84.1866928,88.27938514,91.88055759,
                95.33416644,98.87330272,102.7553798,107.4253971,114.1307069,
                115.0535498,116.0616313,117.1767416,118.4304223,119.8709396,
                121.5779523,123.6981806,126.5538109,131.1412169,132.551708,
                134.3289491,136.7625326,140.7546948,144.5669662,156.2589492,166.8860679],
                [48.614899,79.6328875,85.12588653,89.24148409,92.86224016,
                96.33415833,99.89160385,103.7932662,108.4861933,115.2232423,
                116.150352,117.1630683,118.2832733,119.5426421,120.989644,
                122.7042716,124.8338556,127.7019122,132.3088769,133.7252853,
                135.5099195,137.9535156,141.9618339,145.789233,157.525905,168.1916003],
                [49.31207471,80.54083535,86.06539524,90.20377956,93.84401782,
                97.33415045,100.9098103,104.8309564,109.5466748,116.3152985,
                117.2466535,118.2639807,119.3892541,120.6542817,122.1077349,
                123.8299383,125.9688301,128.8492492,133.4756726,134.897968,
                136.6899578,139.1435155,143.1679079,147.0103582,158.7914944,169.4955745],
                [50.01076179,81.4492545,87.00521409,91.16626852,94.82588913,
                98.3341428,101.9279234,105.8684535,110.6068463,117.4068833,
                118.342462,119.3643766,120.4946926,121.7653501,123.2252218,
                124.9549624,127.1031146,129.9958335,134.6416168,136.0697698,
                137.869078,140.3325473,144.372933,148.2303591,160.0557384,170.7980141],
                [50.71093768,82.35813776,87.94533836,92.12894803,95.80785265,
                99.33413538,102.9459448,106.9057604,111.6667126,118.4980039,
                119.437785,120.4642637,121.5995969,122.8758559,124.3421137,
                126.0793537,128.2367197,131.1416767,135.8067231,137.2407042,
                139.0472942,141.5206256,145.5769251,149.4492527,161.318657,172.098942],
                [51.41258101,83.26747816,88.88576345,93.09181325,96.7899045,
                100.334125,103.9638756,107.9428799,112.7262784,119.5886674,
                120.5326298,121.5636499,122.7039751,123.9858077,125.4584198,
                127.2031216,129.3696555,132.2867897,136.9710038,138.4107841,
                140.2246199,142.7077647,146.7798996,150.6670556,162.5802702,173.3983809],
                [52.1156703,84.17726769,89.82648488,94.0548651,97.77204804,
                101.3341177,104.9817174,108.9798148,113.7855479,120.6788804,
                121.6270035,122.6625424,123.807835,125.0952137,126.5741486,
                128.3262754,130.5019321,133.4311835,138.1344711,139.5800221,
                141.4010683,143.8939785,147.9818717,151.8837839,163.8405971,174.6963526],
                [52.82018479,85.08750198,90.76749664,95.01809908,98.7542797,
                102.3341106,105.9994714,110.0165679,114.8445266,121.7686498,
                122.7209136,123.7609487,124.9111844,126.2040818,127.6893087,
                129.4488243,131.6335592,134.5748687,139.297137,140.7484308,
                142.5766524,145.0792807,149.1828559,153.0994534,165.0996566,175.9928784],
                [53.52610419,85.99817347,91.70879754,95.98151252,99.73659821,
                103.3341037,107.0171388,111.0531418,115.903217,122.8579819,
                123.8143658,124.8588758,126.0140305,127.3124201,128.8039083,
                130.5707769,132.7645461,135.7178556,140.4590131,141.9160221,
                143.7513846,146.2636845,150.3828668,154.3140795,166.3574671,177.2879794],
                [54.23340873,86.90927582,92.65038195,96.94510282,100.7190023,
                104.334097,108.0347211,112.0895391,116.9616241,123.946883,
                124.9073673,125.9563307,127.1163808,128.4202361,129.9179552,
                131.692142,133.8949023,136.8601543,141.620111,143.0828079,
                144.9252774,147.4472029,151.5819183,155.5276771,167.6140464,178.5816757],
                [54.94207905,87.82080288,93.59224581,97.90886745,101.7014907,
                105.3340905,109.052218,113.1257624,118.019752,125.0353594,
                125.9999245,127.05332,128.2182423,129.5275373,131.0314582,
                132.812928,135.0246366,138.0017746,142.7804416,144.2487996,
                146.0983426,148.6298484,152.780024,156.7402609,168.8694119,179.8739874],
                [55.65209626,88.7327486,94.53438513,98.87280394,102.6840598,
                106.3340809,110.0696333,114.1618141,119.0776045,126.1234171,
                127.0920437,128.1498505,129.3196219,130.634331,132.1444244,
                133.933143,136.153758,139.1427261,143.9400158,145.4140085,
                147.270592,149.8116335,153.9771975,157.9518453,170.1235809,181.1649338],
                [56.3634419,89.64510711,95.47679606,99.83690769,103.6667131,
                107.3340744,111.0869668,115.1976966,120.1351857,127.211062,
                128.1837312,129.2459286,130.4205265,131.7406244,133.2568616,
                135.0527951,137.2822749,140.2830183,145.0988443,146.5784456,
                148.442037,150.9925701,155.1734519,159.1624442,171.37657,182.4545339],
                [57.07609793,90.55787143,96.41947479,100.8011804,104.6494473,
                108.3340681,112.1042198,116.2334123,121.1924992,128.2982998,
                129.2749928,130.3415605,131.5209627,132.8464245,134.368777,
                136.1718921,138.4101959,141.4226603,146.2569375,147.7421216,
                149.6126888,152.1726702,156.3687994,160.3720712,172.6283954,183.7428062],
                [57.79004669,91.47103827,97.3624159,101.765618,105.6322611,
                109.334062,113.1213934,117.2689635,122.2495488,129.3851361,
                130.3658345,131.4367525,132.6209368,133.9517379,135.4801778,
                137.2904418,139.5375292,142.561661,147.4143053,148.9050466,
                150.7825583,153.351945,157.5632529,161.5807397,173.8790731,185.029769],
                [58.50527095,92.38460103,98.30561907,102.7302179,106.6151536,
                110.3340561,114.1384886,118.3043525,123.3063381,130.4715765,
                131.456262,132.5315111,133.7204553,135.0565715,136.5910711,
                138.4084514,140.6642827,143.7000292,148.5709578,150.0672315,
                151.9516561,154.530406,158.7568244,162.7884627,175.1286187,186.3154398],
                [59.22175355,93.29855434,99.24907918,103.6949782,107.5981236,
                111.3340504,115.1555065,119.3395814,124.3628706,131.5576263,
                132.546281,133.6258411,134.8195242,136.1609316,137.7014637,
                139.5259285,141.7904644,144.8377735,149.7269046,151.228686,
                153.1199928,155.7080641,159.9495258,163.995253,176.3770473,187.599836],
                [59.93947851,94.21289297,100.1927928,104.6598966,108.5811676,
                112.3340414,116.1724482,120.3746524,125.4191497,132.6432908,
                133.6358969,134.7197488,135.9181497,137.2648247,138.8113624,
                140.64288,142.9160818,145.9749022,150.8821551,152.3894199,
                154.2875784,156.88493,161.1413689,165.2011231,177.6243739,188.8829745],
                [60.65842944,95.12761177,101.1367565,105.6249689,109.5642894,
                113.3340357,117.1893147,121.4095687,126.4751789,133.728575,
                134.7251151,135.8132398,137.0163377,138.3682569,139.9207738,
                141.7593131,144.0411426,147.1114236,152.0367187,153.5494428,
                155.4544229,158.0610143,162.3323649,166.406085,178.870613,190.1648719],
                [61.37859054,96.04270574,102.0809671,106.5901971,110.5474856,
                114.3340301,118.2061069,122.4443302,127.5309614,134.813484,
                135.8139408,136.9063196,138.1140938,139.4712344,141.0297041,
                142.8752345,145.165654,148.247346,153.1906042,154.708764,
                156.6205361,159.2363271,163.522525,167.6101507,180.1157788,191.4455443],
                [62.09994634,96.95816997,103.0254212,107.5555772,111.5307554,
                115.3340247,119.2228259,123.47894,128.5865005,135.8980228,
                136.9023792,137.9989936,139.2114239,140.5737631,142.1381599,
                143.9906509,146.2896232,149.3826764,154.3438207,155.8673926,
                157.7859276,160.4108787,164.7118601,168.813332,181.3598852,192.7250077],
                [62.8224817,97.87399832,103.9701139,108.5211073,112.5140977,
                116.3340194,120.2394725,124.5134001,129.6417993,136.982196,
                137.9904353,139.091267,140.3083333,141.6758487,143.2461471,
                145.1055689,147.4130572,150.5174229,155.4963767,157.0253377,
                158.9506065,161.5846787,165.9003809,170.0156401,182.6029459,194.0032775],
                [63.5461818,98.79018865,104.9150456,109.4867854,113.4975117,
                117.3340143,121.2560478,125.5477124,130.696861,138.0660085,
                139.0781141,140.183145,141.4048275,142.7774971,144.3536719,
                146.2199949,148.5359629,151.6515931,156.6482807,158.1826079,
                160.1145822,162.7577369,167.0880976,171.2170863,183.8449742,195.2803689],
                [64.27103213,99.70673513,105.8602116,110.4526096,114.4809963,
                118.3340059,122.2725525,126.5818787,131.7516884,139.1494648,
                140.1654202,141.2746326,142.5009118,143.8787137,145.46074,
                147.3339346,149.6583471,152.7851941,157.799541,159.3392119,
                161.2778635,163.9300626,168.2750206,172.4176815,185.0859832,196.5562968],
                [64.99701847,100.6236333,106.805609,111.4185781,115.464548,
                119.3340008,123.2889877,127.6159009,132.8062847,140.2325694,
                141.2523585,142.3657348,143.5965913,144.9795041,146.5673574,
                148.4473952,150.7802163,153.9182333,158.9501658,160.4951581,
                162.4404592,165.1016651,169.4611599,173.6174363,186.3259856,197.8310757],
                [65.7241269,101.5408786,107.7512347,112.3846869,116.4481711,
                120.3339958,124.3053542,128.6497808,133.8606527,141.3153267,
                142.3389336,143.4564564,144.6918719,146.0798735,147.6735296,
                149.5603821,151.901577,155.0507176,160.1001629,161.6504546,
                163.6023778,166.2725534,170.6465252,174.8163613,187.5649939,199.1047198],
                [66.45234376,102.458467,108.697086,113.3509384,117.4318623,
                121.3339909,125.3216528,129.6835203,134.9147953,142.3977411,
                143.4251499,144.5468022,145.7867571,147.1798272,148.7792621,
                150.6729014,153.0224354,156.1826539,161.2495401,162.8051096,
                164.7636277,167.4427364,171.831126,176.0144667,188.8030205,200.3772432],
                [67.18165569,103.376394,109.6431602,114.3173289,118.4156207,
                122.3339862,126.3378844,130.7171209,135.9687152,143.4798167,
                144.5110118,145.6367767,146.8812525,148.2793703,149.8845604,
                151.7849586,154.1427978,157.314049,162.3983051,163.9591309,
                165.9242172,168.6122227,173.0149719,177.2117627,190.0400773,201.6486594],
                [67.91204928,104.2946556,110.5894526,115.2838566,119.3994455,
                123.3339816,127.3540511,131.7505845,137.0224151,144.5615578,
                145.5965239,146.7263845,147.9753628,149.3785079,150.9894298,
                152.8965595,155.2626702,158.4449095,163.5464654,165.1125262,
                167.0841542,169.7810208,174.1980719,178.4082589,191.2761761,202.9189819],
                [68.64351226,105.2132462,111.535964,116.25052,120.383336,
                124.3339772,128.3701512,132.7839126,138.0758978,145.6429683,
                146.6816901,147.8156301,149.0690927,150.4772448,152.0938755,
                154.0077097,156.3820586,159.5752418,164.6940282,166.2653032,
                168.2434468,170.949139,175.3804351,179.6039651,192.5113283,204.1882239],
                [69.37603177,106.1321647,112.4826903,117.2173173,121.3672912,
                125.3339691,129.3861867,133.8171069,139.1291658,146.7240523,
                147.7665148,148.9045177,150.1624468,151.5755859,153.1979025,
                155.1184144,157.5009688,160.7050523,165.8410007,167.4174691,
                169.4021021,172.1165855,176.5620703,180.7988907,193.7455454,205.4563981],
                [70.10959544,107.0514059,113.4296287,118.1842469,122.3513075,
                126.3339646,130.4021584,134.8501691,140.1822217,147.8048136,
                148.851002,149.9930518,151.2554296,152.6735359,154.3015159,
                156.2286791,158.6194064,161.8343473,166.98739,168.5690315,
                170.5601287,173.2833683,177.7429862,181.993045,194.9788384,206.7235172],
                [70.84419118,107.9709658,114.3767768,119.151305,123.3353898,
                127.3339603,131.4180669,135.8831006,141.235068,148.8852561,
                149.9351557,151.0812364,152.3480455,153.7710994,155.4047206,
                157.3385089,159.7373771,162.9631327,168.1332029,169.7199972,
                171.7175335,174.4494952,178.9231912,183.1864371,196.2112182,207.9895936],
                [71.5798071,108.8908408,115.3241322,120.1184945,124.3195346,
                128.333956,132.4339132,136.9159031,142.2877061,149.9653827,
                151.0189798,152.1690756,153.4402987,154.868281,156.5075213,
                158.447909,160.8548864,164.0914147,169.2784462,170.8703734,
                172.8743237,175.6149739,180.1026936,184.379076,197.4426953,209.2546394],
                [72.31643155,109.8110272,116.2716924,121.0858118,125.3037412,
                129.3339519,133.4496978,137.948578,143.3401406,151.0451987,
                152.1024782,153.2565736,154.5321935,155.9650851,157.6099228,
                159.5568843,161.9719396,165.219199,170.4231266,172.0201669,
                174.0305066,176.7798118,181.2815017,185.5709702,198.6732803,210.5186665],
                [73.05405309,110.7315214,117.219455,122.0532553,126.2880089,
                130.3339479,134.4654215,138.9811267,144.3923727,152.1247068,
                153.1856545,154.3437342,155.623734,157.0615161,158.7119297,
                160.6654398,163.0885419,166.3464914,171.5672505,173.1693843,
                175.1860892,177.9440165,182.4596234,186.7621284,199.9029834,211.7816865],
                [73.7926605,111.6523198,118.1674158,123.0208236,127.272337,
                131.3339403,135.481085,140.0135508,145.4444048,153.2039107,
                154.2685127,155.4305614,156.7149244,158.1575782,159.8135465,
                161.7735802,164.2046987,167.4732975,172.7108242,174.3180324,
                176.3410783,179.1075952,183.6370665,187.9525591,201.1318146,213.0437109],
                [74.53224279,112.5734176,119.1155763,123.9885152,128.2567247,
                132.3339363,136.496689,141.0458517,146.4962391,154.2828138,
                155.3510561,156.5170588,157.8057685,159.2532764,160.9147777,
                162.8813104,165.3204148,168.5996229,173.8538541,175.4661176,
                177.4954807,180.2705548,184.8138388,189.1422704,202.3597838,214.304751],
                [75.27278913,113.4948141,120.0639323,124.9563289,129.2411684,
                133.3339323,137.5122342,142.0780307,147.547878,155.3614195,
                156.4332885,157.6032303,158.8962704,160.3486134,162.0156276,
                163.9886349,166.4356954,169.725473,174.9963463,176.6136462,
                178.6493029,181.4329025,185.9899478,190.3312705,203.5869006,215.5648177],
                [76.01428891,114.4165046,121.0124817,125.9242632,130.2256732,
                134.3339285,138.5277212,143.1100893,148.5993236,156.4397311,
                157.5152132,158.6890795,159.9864338,161.443594,163.1161004,
                165.0955582,167.5505452,170.8508531,176.1383068,177.7606244,
                179.8025516,182.5946451,187.1654009,191.5195672,204.8131747,216.8239218],
                [76.75673173,115.338486,121.9612224,126.8923142,131.2102358,
                135.3339247,139.5431507,144.1420275,149.650578,157.5177519,
                158.5968338,159.77461,161.0762625,162.538222,164.2162005,
                166.202085,168.6649692,171.9757685,177.2797415,178.9070584,
                180.9552329,183.7557893,188.3402055,192.7071685,206.0386152,218.0820741],
                [77.50010701,116.2607548,122.9101521,127.8604855,132.1948554,
                136.3339211,140.5585232,145.173849,150.7016434,158.5954851,
                159.6781536,160.8598253,162.1657603,163.6325014,165.3159318,
                167.3082194,169.7789719,173.1002242,178.4206562,180.0529542,
                182.1073533,184.9163417,189.5143685,193.8940819,207.2632313,219.3392848],
                [78.24440536,117.1833081,123.859269,128.8287735,133.1795314,
                137.3339176,141.5738395,146.2055541,151.7525218,159.6729339,
                160.7591759,161.9447289,163.2549307,164.7264361,166.4152985,
                168.4139659,170.892558,174.2242254,179.5610567,181.1983177,
                183.2589188,186.0763087,190.6878972,195.080315,208.4870321,220.5955643],
                [78.9896166,118.1061427,124.8085689,129.7971769,134.1642633,
                138.3339103,142.5891002,147.2371439,152.8032154,160.7501014,
                161.839903,163.0293242,164.3437773,165.8200297,167.5143046,
                169.5193286,172.0057321,175.3477769,180.7009484,182.3431546,
                184.4099355,187.2356967,191.8607982,196.2658751,209.710026,221.8509226],
                [79.73573103,119.0292556,125.7580536,130.7656943,135.1490504,
                139.3339067,143.6043058,148.2686196,153.853726,161.8269906,
                162.9203399,164.1136144,165.4323037,166.913286,168.6129538,
                170.6243118,173.1184985,176.4708836,181.8403369,183.4874706,
                185.5604092,188.394512,193.0330785,197.4507694,210.9322223,223.1053696],
                [80.48273914,119.9526423,126.7077193,131.7343247,136.1338888,
                140.3339032,144.6194569,149.2999826,154.9040556,162.9036045,
                164.0004888,165.197603,166.5205132,168.0062085,169.7112501,
                171.7289194,174.2308612,177.5935503,182.9792276,184.6312714,
                186.7103459,189.5527607,194.2047446,198.6350052,212.1536292,224.3589149],
                [81.23063157,120.8763026,127.6575642,132.7030668,137.1187843,
                141.3338998,145.6345542,150.3312341,155.9542063,163.9799461,
                165.0803528,166.281293,167.6084092,169.0988009,170.8091972,
                172.8331555,175.3428254,178.7157817,184.1176258,185.7745623,
                187.8597512,190.7104488,195.3758031,199.8185894,213.374255,225.6115681],
                [81.97939911,121.8002325,128.6075862,133.6719168,138.1037332,
                142.3338965,146.6495981,151.3623751,157.0041798,165.0560182,
                166.159935,167.3646877,168.6959951,170.1910667,171.9067988,
                173.937024,176.454395,179.8375823,185.2555366,186.9173488,
                189.0086307,191.8675822,196.5462603,201.0015291,214.5941081,226.8633385],
                [82.72903274,122.7244289,129.5577836,134.6408786,139.088735,
                143.3338933,147.6645893,152.3934069,158.0539781,166.1318237,
                167.2392382,168.4477901,169.7832741,171.2830092,173.0040586,
                175.0405288,177.5655739,180.9589568,186.3929652,188.0596362,
                190.1569899,193.0241668,197.7161226,202.1838306,215.8131963,228.1142354],
                [83.47952358,123.6488892,130.5081545,135.6099485,140.0737892,
                144.3338902,148.6795266,153.4243306,159.1036029,167.2073653,
                168.3182655,169.5306034,170.8702495,172.3746318,174.1009801,
                176.1436736,178.6763664,182.0799095,187.5299167,189.2014296,
                191.3048343,194.1802083,198.8853962,203.3655007,217.0315278,229.3642678],
                [84.23086289,124.5736106,131.4586971,136.5791256,141.0588951,
                145.3338831,149.6944138,154.4551474,160.1530561,168.2826458,
                169.3970196,170.6131304,171.9569243,173.4659379,175.1975667,
                177.2464623,179.7867763,183.2004448,188.6663959,190.3427343,
                192.4521692,195.3357122,200.0540871,204.5465459,218.2491102,230.6134445],
                [84.98304208,125.4985904,132.4094096,137.5484086,142.0440523,
                146.3338799,150.7092499,155.4858583,161.2023395,169.3576679,
                170.4755033,171.6953743,173.0433016,174.5569307,176.2938221,
                178.3488984,180.8968077,184.320567,189.8024077,191.4835551,
                193.5989997,196.4906841,201.2222014,205.7269727,219.4659512,231.8617743],
                [85.73605272,126.423826,133.3602882,138.5177965,143.0292571,
                147.3338768,151.7240354,156.5164644,162.2514547,170.4324342,
                171.5537196,172.7773378,174.1293846,175.6476134,177.3897494,
                179.4509856,182.0064644,185.4402803,190.9379569,192.6238971,
                194.7453312,197.6451294,202.3897448,206.9067873,220.6820584,233.1092658],
                [86.48988651,127.3493147,134.3113353,139.4872882,144.014515,
                148.3338738,152.7387708,157.5469668,163.3004036,171.5069473,
                172.6316709,173.8590237,175.2151761,176.7379892,178.4853522,
                180.5527273,183.1157502,186.5595888,192.0730482,193.763765,
                195.8911685,198.7990535,203.5567232,208.085996,221.8974392,234.3559274],
                [87.24453528,128.2750525,135.2625471,140.4568827,144.9998226,
                149.3338708,153.7534566,158.5773665,164.3491877,172.5812097,
                173.7093601,174.940435,176.3006792,177.8280612,179.5806343,
                181.6541272,184.2246689,187.6784968,193.2076861,194.9031637,
                197.0365168,199.9524615,204.7231423,209.2646048,223.1121008,235.6017674],
                [87.99999062,129.2010399,136.2139222,141.4265763,145.9851795,
                150.333868,154.7680934,159.6076646,165.3978087,173.6552241,
                174.7867898,176.0215735,177.3858967,178.9178323,180.6755976,
                182.7551886,185.3332241,188.797008,194.3418752,196.0420979,
                198.1813808,201.1053587,205.8890075,210.4426198,224.3260504,236.8467939],
                [88.75624537,130.1272729,137.1654587,142.3963731,146.970585,
                151.3338652,155.7826815,160.6378621,166.4462682,174.7289928,
                175.8639625,177.1024436,178.4708315,180.0073057,181.7702459,
                183.8559148,186.4414194,189.9151266,195.47562,197.1805721,
                199.3257654,202.2577501,207.0543246,211.6200467,225.5392951,238.0910151],
                [89.5132914,131.0537491,138.1171552,143.3662696,147.9560388,
                152.3338584,156.7972215,161.6679599,167.4945679,175.8025183,
                176.9408808,178.1830471,179.5554863,181.0964842,182.8645824,
                184.9563093,187.5492584,191.0328563,196.6089249,198.318591,
                200.4696753,203.4096408,208.2190987,212.7968915,226.7518418,239.3344388],
                [90.27112106,131.9804661,139.0690101,144.3362649,148.9415404,
                153.3338555,157.8117138,162.6979591,168.5427093,176.875803,
                178.0175472,179.2633866,180.6398639,182.1853708,183.95861,
                186.0563752,188.6567446,192.1502011,197.7417941,199.4561589,
                201.6131152,204.5610351,209.3833353,213.9731597,227.9636972,240.5770727],
                [91.02972684,132.9074216,140.0210217,145.3063579,149.9270893,
                154.3338527,158.826159,163.7278606,169.590694,177.9488492,
                179.0939641,180.3434648,181.7239671,183.2739682,185.0523319,
                187.1561157,189.7638814,193.2671645,198.874232,200.5932803,
                202.7560896,205.7119387,210.5470396,215.148857,229.1748681,241.8189245],
                [91.78910132,133.8346133,140.9731865,146.2765478,150.9126816,
                155.33385,159.8405574,164.7576653,170.6385234,179.0216592,
                180.170134,181.4232841,182.8077984,184.3622793,186.145751,
                188.2555341,190.8706722,194.3837504,200.0062426,201.7299594,
                203.8986031,206.8623561,211.7102166,216.3239889,230.3853611,243.0600018],
                [92.54923723,134.7620388,141.9255071,147.2468336,151.8983234,
                156.3338474,160.8549096,165.7873743,171.6861992,180.0942354,
                181.2460592,182.502847,183.8913605,185.4503069,187.2388701,
                189.3546334,191.9771204,195.4999624,201.1378301,202.8662006,
                205.04066,208.0122918,212.8728716,217.4985609,231.5951827,244.3003119],
                [93.31012737,135.689696,142.87798,148.2172144,152.8840112,
                157.3338448,161.8692159,166.8169883,172.7337226,181.1665799,
                182.3217421,183.5821561,184.974656,186.5380537,188.3316923,
                190.4534166,193.0832292,196.615804,202.2689986,204.002008,
                206.1822648,209.1617505,214.0350094,218.6725782,232.8043392,245.539862],
                [94.0717647,136.617581,143.8306037,149.1876865,153.8697445,
                158.3338424,162.8834769,167.8465084,173.7810966,182.2386949,
                183.397185,184.6612137,186.0576874,187.6255222,189.4242202,
                191.5518868,194.1890018,197.7312788,203.399752,205.1373857,
                207.3234217,210.3107366,215.196635,219.8460462,234.012837,246.7786595],
                [94.83414227,137.5456948,144.7833769,150.1582543,154.8555229,
                159.3338358,163.8976929,168.8759353,174.8283199,183.3105827,
                184.4723902,185.7400221,187.1404572,188.7127153,190.5164567,
                192.6500469,195.2944415,198.8463902,204.5300942,206.2723379,
                208.4641349,211.4592547,216.3577532,221.01897,235.2206822,248.0167113],
                [95.59725322,138.4740337,145.7362982,151.1289145,155.8413459,
                160.3338332,164.9118643,169.9052701,175.8753953,184.3822453,
                185.5473598,186.8185838,188.2229679,189.7996354,191.6084045,
                193.7478999,196.3995513,199.9611417,205.6600291,207.4068684,
                209.6044087,212.6073091,217.5183687,222.1913546,236.427881,249.2540244],
                [96.36109082,139.4025956,146.6893661,152.0996661,156.8272132,
                161.3338307,165.9259916,170.934515,176.9223241,185.4536849,
                186.6220962,187.8969011,189.3052219,190.8862851,192.7000664,
                194.8454486,197.5043344,201.0755365,206.7895605,208.5409812,
                210.744247,213.7549041,218.6784862,223.3632052,237.6344392,250.4906056],
                [97.12564843,140.3313785,147.6425792,153.0705083,157.8131209,
                162.3338283,166.9400752,171.9636679,177.9691077,186.5249036,
                187.6966014,188.9749761,190.3872215,191.9726669,193.7914449,
                195.9426958,198.6087937,202.189578,207.9186921,209.6746801,
                211.883654,214.9020441,219.8381103,224.5345266,238.8403629,251.7264616],
                [97.89091951,141.2603803,148.5959342,154.0414403,158.7990751,
                163.3338259,167.9541154,172.9927312,179.0157474,187.5959033,
                188.7708776,190.0528112,191.4689691,193.0587832,194.8825426,
                197.0396443,199.7129322,203.3032695,209.0474276,210.8079691,
                213.0226336,216.0487331,220.9972456,225.7053236,240.0456579,252.9615991],
                [98.65689721,142.1895991,149.5494338,155.0124612,159.7850724,
                164.3338236,168.9681128,174.0217056,180.0622445,188.666686,
                189.8449268,191.1304085,192.5504662,194.1446366,195.9733622,
                198.1362968,200.8167529,204.4166141,210.1757707,211.9408518,
                214.1611897,217.1949753,222.1558966,226.8756012,241.2503298,254.1960247],
                [99.42357599,143.1190328,150.5030747,155.9835703,160.7711123,
                165.3338214,169.9820676,175.050592,181.1086003,189.7372539,
                190.9187513,192.2077703,193.6317167,195.2302294,197.0639062,
                199.232656,201.9202587,205.5296151,211.3037248,213.0733319,
                215.2993261,218.3407747,223.3140675,228.045364,242.4543844,255.4297446],
                [100.1909492,144.0486794,151.4568557,156.9547638,161.7571944,
                166.3338192,170.9959821,176.0793911,182.1548162,190.8076087,
                191.9923529,193.2848985,194.712722,196.315564,198.154177,
                200.3287247,203.0234523,206.6422754,212.4312935,214.2054131,
                216.4370467,219.4861354,224.4717628,229.2146165,243.6578271,256.6627654],
                [100.9590107,144.9785355,152.4107755,157.9260465,162.7433184,
                167.3338127,172.0098531,177.1081038,183.2008933,191.8777525,
                193.0657337,194.3617954,195.7934844,197.4006426,199.2441772,
                201.4245053,204.1263368,207.7545983,213.5584809,215.3370989,
                217.5743552,220.6310613,225.6289867,230.3833636,244.8606635,257.8950932],
                [101.7277544,145.9086024,153.3648328,158.897415,163.7294839,
                168.3338104,173.0236828,178.1367308,184.246833,192.947687,
                194.1388957,195.438463,196.8740059,198.4854676,200.3339091,
                202.5200005,205.2289147,208.8665866,214.6852892,216.4683929,
                218.7112552,221.7755563,226.7857434,231.5516095,246.0628989,259.1267342],
                [102.4971744,146.8388766,154.3190264,159.8688685,164.7156904,
                169.3338082,174.0374714,179.1652729,185.2926364,194.0174142,
                195.2118407,196.5149034,197.9542889,199.5700412,201.4233752,
                203.6152129,206.3311888,209.9782434,215.8117223,217.5992985,
                219.8477504,222.9196241,227.942037,232.7193588,247.2645386,260.3576944],
                [103.2672647,147.7693562,155.2733552,160.8404063,165.7019343,
                170.3338061,175.0512194,180.1937308,186.3383047,195.086937,
                196.2845709,197.5911185,199.0343352,200.6543657,202.5125777,
                204.7101448,207.4331618,211.0895716,216.9377836,218.7298193,
                220.9838444,224.0632686,229.0978718,233.8866158,248.465588,261.5879798],
                [104.0380195,148.7000396,156.2278178,161.8120275,166.6882218,
                171.333804,176.064927,181.2221053,187.3838392,196.1562551,
                197.3570879,198.6671103,200.1141472,201.7384432,203.601519,
                205.8047988,208.5348363,212.200574,218.0634762,219.8599584,
                222.1195406,225.2064935,230.2532515,235.0533848,249.666052,262.8175963],
                [104.8094332,149.6309249,157.1824109,162.7837316,167.6745493,
                172.333802,177.0785947,182.250397,188.429241,197.2253712,
                198.4293938,199.7428809,201.1937267,202.822276,204.6902014,
                206.8991781,209.636215,213.3112536,219.1888036,220.9897194,
                223.2548425,226.3493024,231.4081803,236.2196701,250.8659359,264.0465496],
                [105.5815,150.5620103,158.1371378,163.7555177,168.6609165,
                173.3338001,178.0922229,183.2786068,189.4745112,198.2942872,
                199.5014902,200.818432,202.2730758,203.905866,205.7786271,
                207.9932835,210.7373003,214.4216131,220.3137687,222.1191054,
                224.3897537,227.4916989,232.562662,237.3854759,252.0652447,265.2748456],
                [106.3542145,151.4932941,159.0919952,164.7273822,169.647323,
                174.3337982,179.1058117,184.3067353,190.519651,199.3630047,
                200.5733791,201.8937656,203.3521965,204.9892155,206.8667984,
                209.0871182,211.838095,215.5316553,221.4383747,223.2481197,
                225.5242773,228.6336867,233.7167005,238.5508063,253.2639832,266.5024898],
                [107.1275711,152.4247747,160.0469818,165.6993302,170.6337686,
                175.3337917,180.1193617,185.3347833,191.5646615,200.4315256,
                201.6450622,202.9688837,204.4310907,206.0723264,207.9547173,
                210.1806844,212.9386013,216.6413829,222.5626247,224.3767655,
                226.6584168,229.775269,234.8702995,239.7156654,254.4621563,267.7294878],
                [107.9015644,153.3564503,161.0020967,166.671358,171.6202528,
                176.3337897,181.1328731,186.3627513,192.6095439,201.4998513,
                202.7165413,204.0437879,205.5097605,207.1552009,209.0423862,
                211.2739845,214.0388219,217.7507986,223.6865218,225.5050459,
                227.7921754,230.9164496,236.0234629,240.8800571,255.6597689,268.9558451],
                [108.6761891,154.2883176,161.9573388,167.6434652,172.6067754,
                177.3337878,182.1463462,187.3906402,193.6542991,202.5679837,
                203.787818,205.1184801,206.5882076,208.2378409,210.1298071,
                212.3670208,215.1387591,218.859905,224.8100688,226.632964,
                228.9255563,232.0572316,237.1761943,242.0439855,256.8568257,270.1815671],
                [109.45144,155.2203782,162.9127069,168.6156509,173.5933324,
                178.3337859,183.1597813,188.4184505,194.6989283,203.6359244,
                204.8588942,206.1929622,207.666434,209.3202484,211.216982,
                213.4597953,216.2384152,219.9687047,225.9332689,227.7605229,
                230.0585628,233.1976186,238.3284973,243.2074543,258.0533313,271.406659],
                [110.2273113,156.1526291,163.8682,169.5879147,174.5799305,
                179.333784,184.1731789,189.4461829,195.7434326,204.7036749,
                205.9297715,207.2672358,208.7444416,210.4024254,212.3039131,
                214.5523104,217.3377928,221.0772004,227.0561248,228.8877255,
                231.1911979,234.3376137,239.4803755,244.3704676,259.2492903,272.6311261],
                [111.003799,157.0850684,164.8238171,170.5602557,175.5665659,
                180.3337823,185.1865392,190.4738381,196.787813,205.7712368,
                207.0004515,208.3413028,209.8222321,211.4843737,213.3906024,
                215.6445683,218.436894,222.1853938,228.1786394,230.0145749,
                232.3234648,235.4772203,240.6318324,245.533029,260.4447074,273.8549737],
                [111.7808974,158.0176948,165.7795548,171.5326733,176.5532385,
                181.3337806,186.1998624,191.5014167,197.8320704,206.8386117,
                208.0709359,209.4151648,210.8998073,212.5660953,214.4770519,
                216.7365709,219.5357213,223.2932889,229.3008157,231.141074,
                233.4553665,236.6164416,241.7828716,246.6951423,261.6395868,275.0782068],
                [112.5586017,158.9505066,166.7354166,172.505167,177.5399478,
                182.3337789,187.213149,192.5289194,198.876206,207.9058012,
                209.1412263,210.4888235,211.9771691,213.6475921,215.5632636,
                217.8283205,220.6342767,224.4008873,230.4226563,232.2672256,
                234.586906,237.7552808,242.9334965,247.8568112,262.8339332,276.3008305],
                [113.3369068,159.8835023,167.6913994,173.4777331,178.5266935,
                183.3337725,188.2263992,193.5563467,199.9202207,208.9728068,
                210.2113243,211.5622807,213.0543192,214.7288659,216.6492394,
                218.9198192,221.7325627,225.5081916,231.544164,233.3930327,
                235.7180863,238.8937409,244.0837104,249.0180392,264.0277508,277.5228496],
                [114.1158081,160.8166805,168.6475021,174.4503768,179.5134754,
                184.3337707,189.2396133,194.5836993,200.9641156,210.0396299,
                211.2812326,212.6355378,214.1312593,215.8099175,217.7349813,
                220.0110689,222.8305813,226.6152042,232.6653416,234.5184979,
                236.8489103,240.0318251,245.2335168,250.1788301,265.2210438,278.7442692],
                [114.8953006,161.7500396,169.6037237,175.4230946,180.5002932,
                185.333769,190.2527916,195.6109778,202.0078915,211.1062721,
                212.3509505,213.7085967,215.2079911,216.8907507,218.8204911,
                221.1020717,223.9283348,227.7219273,233.7861916,235.6436241,
                237.9793808,241.1695365,246.3829189,251.3391874,266.4138167,279.965094],
                [115.6753797,162.6835782,170.5600633,176.3958861,181.4871466,
                186.3337673,191.2659345,196.6381828,203.0515494,212.1727348,
                213.4204807,214.7814588,216.2845163,217.9713662,219.9057707,
                222.1928296,225.0258252,228.8283634,234.9067168,236.7684141,
                239.1095009,242.3068779,247.5319199,252.4991144,267.6060734,281.1853289],
                [116.4560407,163.6172929,171.51652,177.3687505,182.4740314,
                187.3337657,192.2790421,197.6653149,204.0950904,213.2390194,
                214.4898246,215.8541257,217.3608365,219.0517659,220.9908221,
                223.2833445,226.1230547,229.9345148,236.0269198,237.8928705,
                240.2392732,243.4438525,248.6805231,253.6586148,268.7978182,282.4049785],
                [117.2372791,164.5511861,172.4730928,178.3416873,183.4609549,
                188.3337642,193.2921147,198.6923746,205.1385153,214.3051275,
                215.5589838,216.9265991,218.4369535,220.1319515,222.0756469,
                224.3736185,227.2200255,231.0403837,237.1468031,239.0169959,
                241.3687005,244.580463,249.8287317,254.8176918,269.9890552,283.6240475],
                [118.0190903,165.4852544,173.4297807,179.3146959,184.4479131,
                189.3337627,194.3051527,199.7193626,206.1818251,215.3710603,
                216.6279597,217.9988804,219.5128688,221.2119246,223.160247,
                225.4636533,228.3167394,232.1459723,238.2663692,240.140793,
                242.4977855,245.7167125,250.9765488,255.9763488,271.1797883,284.8425404],
                [118.8014699,166.4194966,174.386583,180.2877758,185.4349057,
                190.3337612,195.3181563,200.7462777,207.2250206,216.4368194,
                217.6967538,219.0709712,220.5885839,222.2916869,224.2446242,
                226.553451,229.4131986,233.251283,239.3856207,241.2642645,
                243.6265311,246.8526037,252.1239774,257.1345892,272.3700214,286.0604618],
                [119.5844133,167.3539112,175.343496,181.2609264,186.4219326,
                191.3337549,196.3311237,201.7731238,208.2681028,217.5024059,
                218.7653675,220.142873,221.6641006,223.3712402,225.3287803,
                227.6430133,230.5094051,234.3563178,240.5045601,242.3874134,
                244.7549398,247.9881394,253.2710207,258.2924161,273.5597586,287.2778161],
                [120.3679163,168.288497,176.300524,182.2341439,187.4089933,
                192.3337533,197.3440593,202.7998998,209.3110726,218.5678214,
                219.8338021,221.2145872,222.7394203,224.450586,226.4127169,
                228.7323421,231.6053608,235.4610789,241.6231897,243.5102411,
                245.8830142,249.1233225,254.4176815,259.449833,274.7490037,288.4946078],
                [121.1519746,169.2232525,177.2576636,183.2074341,188.3960878,
                193.3337518,198.3569613,203.8266063,210.3539307,219.633067,
                220.9020592,222.2861153,223.8145446,225.5297259,227.4964357,
                229.8214393,232.7010677,236.5655685,242.741512,244.6327506,
                247.010757,250.2581557,255.563963,260.6068428,275.9377605,289.7108412],
                [121.9365839,170.1581765,178.2149139,184.1807932,189.3832156,
                194.3337503,199.3698299,204.8532438,211.3966767,220.6981442,
                221.9701401,223.3574587,224.8894749,226.6086616,228.5799386,
                230.9103065,233.7965278,237.6697887,243.8595294,245.7549445,
                248.1381708,251.3926417,256.709868,261.7634487,277.1260327,290.9265206],
                [122.7217395,171.0932676,179.1722741,185.1542209,190.3703727,
                195.3337488,200.3826655,205.8798128,212.4393143,221.7630543,
                223.0380461,224.428619,225.9642129,227.6873945,229.663227,
                231.9989457,234.8917428,238.7737414,244.9772441,246.8768251,
                249.265258,252.5267831,257.8553994,262.9196539,278.313824,292.1416502],
                [123.5074383,172.0285246,180.1297434,186.1277166,191.3575663,
                196.3337474,201.3954682,206.9063139,213.4818429,222.8277985,
                224.1057785,225.4995974,227.0387599,228.7659264,230.7463027,
                233.0873585,235.9867148,239.8774289,246.0946586,247.9983949,
                250.3920211,253.6605825,259.00056,264.0754615,279.5011382,293.3562342],
                [124.2936757,172.9639443,181.0873208,187.1012798,192.3447925,
                197.3337461,202.4082384,207.9327476,214.5242633,223.892378,
                225.1733388,226.5703953,228.1131174,229.8442587,231.8291673,
                234.1755466,237.0814455,240.980853,247.2117751,249.1196562,
                251.5184626,254.7940427,260.1453528,265.2308743,280.6879789,294.5702768],
                [125.0804476,173.8995292,182.0450057,188.0749099,193.3320512,
                198.3337449,203.4209763,208.9591144,215.5665763,224.9567942,
                226.2407282,227.6410142,229.1872868,230.9223928,232.9118224,
                235.2635118,238.1759368,242.0840158,248.3285959,250.2406114,
                252.6445851,255.9271661,261.2897804,266.3858955,281.8743495,295.783782],
                [125.8677501,174.8352763,183.0027972,189.0486065,194.3193419,
                199.3337386,204.4336821,209.9854147,216.6087828,226.0210483,
                227.3079479,228.7114554,230.2612696,232.0003304,233.9942695,
                236.3512558,239.2701904,243.1869193,249.4451232,251.3612627,
                253.7703907,257.0599554,262.4338457,267.540528,283.0602538,296.9967539]]

        def __init__(self):
                pass

class lmtbl_ttest:
        seeds=[0.9999,0.999,0.99,0.9,0.8,0.7,0.6,0.5,0.4,0.3,
        0.2,0.1,0.09,0.08,0.07,0.06,0.05,0.04,0.03,0.02,0.01,
        0.008,0.006,0.004,0.002,0.001,0.0005,0.0001,0.00001]
        data=[0.000125131,0.001253336,0.012533475,0.125661379,0.25334717,
        0.385320577,0.52440068,0.674489996,0.841621593,1.036433927,
        1.281552411,1.64485515,1.695399352,1.750687769,1.81191252,
        1.880795634,1.959966301,2.053751528,2.170093411,2.326351582,
        2.575834209,2.652075109,2.747787245,2.878168412,3.090240447,
        3.290536458,3.480767809,3.890607581,4.417196062]
        data1=[
                [0.00015708,0.001570798,0.015709255,0.15838444,0.324919696,
                0.50952545,0.726542528,1,1.37638192,1.962610505,3.077683537,
                6.313751514,7.026366228,7.915815087,9.057886685,10.57889499,
                12.70620473,15.89454484,21.20494879,31.82051595,63.65674115,
                127.3213364,636.6192487,6366.19767,63661.97722],
                [0.000141421,0.001414214,0.014142843,0.142133811,0.288675135,
                0.44474959,0.6172134,0.816496581,1.060660172,1.386206559,
                1.885618083,2.91998558,3.103976694,3.319764048,3.578246638,
                3.89642536,4.30265273,4.848732213,5.642778353,6.964556734,
                9.9248432,14.08904727,31.59905458,99.99249984,316.2253943],
                [0.000136035,0.00136035,0.013604055,0.136598199,0.276670662,
                0.424201622,0.584389727,0.764892328,0.978472312,1.249778105,
                1.637744352,2.353363435,2.470806799,2.605426823,2.762598961,
                2.950510469,3.182446305,3.48190876,3.896045934,4.540702858,
                5.840909309,7.453318505,12.92397864,28.00013001,60.39682404],
                [0.000133333,0.001333334,0.013333827,0.133830367,0.270722295,
                0.41416326,0.568649063,0.740697084,0.940964577,1.189566853,
                1.533206273,2.131846782,2.226099555,2.332872559,2.455891992,
                2.600761995,2.776445105,2.998527873,3.297629727,3.746947388,
                4.604094871,5.597568367,8.610301581,15.54410058,27.77164766],
                [0.000131715,0.001317153,0.013171985,0.132175175,0.267180866,
                0.408228733,0.559429645,0.726686844,0.91954378,1.155767344,
                1.475884037,2.015048372,2.097836667,2.190958257,2.297392326,
                2.421584707,2.570581835,2.756508522,3.002874974,3.364929997,
                4.032142983,4.773340605,6.868826626,11.17771007,17.89686614],
                [0.000130639,0.001306395,0.013064379,0.131075653,0.264834533,
                0.404313361,0.553380924,0.717558197,0.905703285,1.134156931,
                1.439755747,1.943180274,2.019200795,2.104306121,2.201058928,
                2.313263299,2.446911846,2.612241845,2.828927862,3.142668403,
                3.70742802,4.316827103,5.958816179,9.082346327,13.55531642],
                [0.000129873,0.001298731,0.012987719,0.130292797,0.263166861,
                0.401538232,0.549109658,0.711141778,0.896029644,1.119159128,
                1.414923928,1.894578604,1.966152951,2.046011102,2.136452899,
                2.240879287,2.364624251,2.516752418,2.714573009,2.997951566,
                3.499483297,4.029337177,5.407882521,7.884584262,11.2148498],
                [0.0001293,0.001292996,0.012930358,0.129707272,0.261921097,
                0.399469297,0.545933764,0.706386613,0.888889518,1.108145445,
                1.39681531,1.859548033,1.92798552,2.00415154,2.090166012,
                2.189154805,2.306004133,2.448984989,2.633814372,2.896459446,
                3.355387331,3.832518684,5.041305433,7.120003882,9.782531715],
                [0.000128854,0.001288544,0.012885835,0.129252932,0.260955336,
                0.397867765,0.543480241,0.702722147,0.88340386,1.099716197,
                1.383028739,1.833112923,1.899221814,1.972652649,2.055394858,
                2.150375268,2.262157158,2.398440983,2.573803977,2.821437921,
                3.249835541,3.689662392,4.780912586,6.593682584,8.827483617],
                [0.000128499,0.001284989,0.012850279,0.128890189,0.26018483,
                0.396591494,0.541528039,0.699812061,0.879057829,1.093058074,
                1.372183641,1.812461102,1.876774376,1.948099462,2.028327014,
                2.120233532,2.228138842,2.35931462,2.527484243,2.763769458,
                3.169272672,3.581406202,4.586893858,6.211050891,8.150286559],
                [0.000128208,0.001282085,0.012821233,0.128593911,0.25955586,
                0.395550624,0.539937879,0.697445328,0.875529978,1.08766638,
                1.363430318,1.795884814,1.858771961,1.928426823,2.006662752,
                2.096138836,2.200985159,2.328139826,2.490663931,2.718079183,
                3.105806514,3.496614173,4.436979338,5.921194162,7.647487116],
                [0.000127967,0.001279669,0.01279706,0.128347381,0.259032746,
                0.394685581,0.538617668,0.695482866,0.872609292,1.083211421,
                1.356217334,1.782287548,1.844015113,1.912313304,1.988933504,
                2.076440623,2.178812827,2.302721671,2.460700162,2.680997992,
                3.054539586,3.428444242,4.317791282,5.694465793,7.260633919],
                [0.000127763,0.001277626,0.012776629,0.12813905,0.258590858,
                0.393955314,0.53750409,0.693829304,0.870151534,1.079468737,
                1.350171289,1.770933383,1.831699761,1.898874471,1.974157952,
                2.060038063,2.160368652,2.281603562,2.435845205,2.650308836,
                3.012275833,3.37246794,4.220831726,5.512515049,6.954441172],
                [0.000127588,0.001275877,0.012759136,0.127960685,0.258212654,
                0.393330622,0.53655218,0.69241707,0.868054782,1.076280245,
                1.345030375,1.761310115,1.821266949,1.887496138,1.961655674,
                2.046169055,2.144786681,2.263781276,2.414897716,2.624494064,
                2.976842734,3.325695816,4.140454113,5.363413041,6.706452292],
                [0.000127436,0.001274362,0.012743989,0.127806258,0.257885301,
                0.392790166,0.535729133,0.691196949,0.866244973,1.073531396,
                1.340605608,1.753050325,1.812316083,1.877738655,1.950940065,
                2.034289392,2.131449536,2.248540287,2.397005036,2.60248029,
                2.946712883,3.286038568,4.072765196,5.239088212,6.50174629],
                [0.000127304,0.001273038,0.012730746,0.127671257,0.257599195,
                0.392317993,0.535010453,0.690132254,0.864667002,1.071137163,
                1.336757167,1.745883669,1.80455261,1.869279026,1.941654062,
                2.024000177,2.119905285,2.235358425,2.381545371,2.583487179,
                2.920781621,3.251992871,4.014996327,5.133893517,6.33004791],
                [0.000127187,0.001271871,0.01271907,0.127552234,0.257347006,
                0.391901938,0.53437748,0.689195075,0.863279018,1.069033111,
                1.33337939,1.739606716,1.797755158,1.861874664,1.933529711,
                2.01500233,2.109815559,2.223845299,2.368054758,2.566933975,
                2.898230518,3.222449906,3.965126272,5.043764976,6.184064371],
                [0.000127083,0.001270834,0.012708698,0.127446513,0.257123043,
                0.39153256,0.533815751,0.688363807,0.862048668,1.067169516,
                1.330390944,1.734063592,1.791754065,1.855339851,1.926362038,
                2.007067314,2.100922037,2.21370324,2.356180003,2.552379618,
                2.878440471,3.196574222,3.921645825,4.965706285,6.058484531],
                [0.000126991,0.001269907,0.012699423,0.127351983,0.25692282,
                0.391202426,0.533313882,0.68762146,0.86095055,1.065507399,
                1.327728209,1.729132792,1.786417229,1.849530034,1.919991578,
                2.000017464,2.09302405,2.204701336,2.345647528,2.539483189,
                2.860934604,3.17372453,3.883405852,4.897461588,5.949353353],
                [0.000126907,0.001269073,0.012691081,0.127266957,0.256742754,
                0.390905598,0.532862792,0.686954497,0.85996444,1.064015771,
                1.325340707,1.724718218,1.781640205,1.844330934,1.914292421,
                1.993712596,2.085963441,2.196657727,2.336242153,2.527977001,
                2.845339707,3.153400532,3.849516274,4.837301152,5.853667922],
                [0.000126832,0.001268318,0.012683536,0.12719007,0.256579948,
                0.39063728,0.532455147,0.686351991,0.859074035,1.062669688,
                1.323187874,1.720742871,1.777339344,1.839651131,1.909163829,
                1.988040625,2.079613837,2.189427267,2.327792316,2.517648014,
                2.831359554,3.135206245,3.819277164,4.783877115,5.76910899],
                [0.000126763,0.001267633,0.012676681,0.127120209,0.256432034,
                0.390393553,0.532084961,0.685805032,0.858266052,1.061448843,
                1.321236742,1.717144335,1.773446864,1.835416563,1.904524282,
                1.982910864,2.073873058,2.182892646,2.320159567,2.50832455,
                2.818756056,3.118824206,3.792130671,4.736124059,5.693858041],
                [0.000126701,0.001267007,0.012670425,0.127056451,0.25629706,
                0.390171186,0.531747299,0.685306278,0.857529554,1.06033654,
                1.31946024,1.713871517,1.769907261,1.831566603,1.900306998,
                1.978249153,2.068657599,2.176958106,2.313230947,2.499866736,
                2.807335678,3.103996962,3.767626803,4.693188999,5.626469823],
                [0.000126643,0.001266434,0.012664692,0.126998032,0.256173398,
                0.38996749,0.531438056,0.684849627,0.856855458,1.059318921,
                1.317835934,1.710882067,1.766674647,1.828051152,1.896456886,
                1.97399426,2.063898547,2.17154467,2.306913389,2.492159469,
                2.796939498,3.090513547,3.745398618,4.654381147,5.565781762],
                [0.000126591,0.001265907,0.012659419,0.126944307,0.256059685,
                0.389780209,0.53115379,0.684429965,0.856236158,1.058384393,
                1.316345073,1.708140745,1.76371076,1.824828445,1.892928014,
                1.97009521,2.059538536,2.166586627,2.301129524,2.48510717,
                2.787435805,3.078199459,3.725143948,4.619135235,5.510848411],
                [0.000126542,0.001265421,0.012654554,0.126894733,0.255954766,
                0.389607434,0.530891592,0.684042973,0.855665233,1.057523179,
                1.314971864,1.705617901,1.760983475,1.82186339,1.889681804,
                1.966509131,2.055529418,2.162028865,2.295814506,2.478629817,
                2.778714523,3.066909114,3.706611742,4.586984348,5.460893266],
                [0.000126497,0.00126497,0.012650051,0.126848847,0.255857659,
                0.389447545,0.530648989,0.683684979,0.855137231,1.056726981,
                1.313702913,1.703288423,1.758465498,1.819126287,1.886685613,
                1.963199826,2.051830493,2.157824814,2.290913604,2.472659904,
                2.770682946,3.056520106,3.689591711,4.557539482,5.415272805],
                [0.000126455,0.001264553,0.012645871,0.126806251,0.255767523,
                0.389299152,0.530423865,0.683352843,0.854647486,1.055988703,
                1.312526782,1.701130908,1.756133656,1.816591829,1.883911639,
                1.960136457,2.048407115,2.153934856,2.286380225,2.467140089,
                2.763262442,3.046928772,3.673906398,4.530473997,5.373449317],
                [0.000126416,0.001264163,0.01264198,0.126766606,0.255683635,
                0.389161057,0.530214396,0.683043861,0.854191986,1.055302249,
                1.311433647,1.699126996,1.753968051,1.814238321,1.881336051,
                1.957292599,2.045229611,2.150325074,2.282174556,2.46202135,
                2.756385902,3.038046742,3.659405017,4.505511631,5.334970094],
                [0.00012638,0.0012638,0.012638349,0.126729613,0.255605365,
                0.389032226,0.530019004,0.682755693,0.853767262,1.054662347,
                1.310415025,1.697260851,1.751951521,1.812047096,1.878938309,
                1.954645481,2.042272449,2.146966263,2.278262321,2.457261531,
                2.749995652,3.02979822,3.645958632,4.482417175,5.299451342],
                [0.000126346,0.001263461,0.012634953,0.126695016,0.255532168,
                0.388911756,0.529836317,0.682486306,0.853370296,1.054064419,
                1.309463549,1.695518742,1.750069191,1.810001869,1.876700611,
                1.952175377,2.039513438,2.14383314,2.274613858,2.45282418,
                2.744041917,3.02211783,3.633456346,4.460989141,5.26656559],
                [0.000126314,0.001263143,0.01263177,0.12666259,0.255463567,
                0.388798859,0.529665134,0.682233921,0.852998454,1.053504465,
                1.308572793,1.693888703,1.748308098,1.808088549,1.874607455,
                1.949865111,2.036933334,2.140903704,2.271203372,2.448677619,
                2.73848148,3.014948884,3.621802256,4.441053938,5.236031765],
                [0.000126284,0.001262844,0.012628781,0.126632135,0.255399142,
                0.388692843,0.529504403,0.681996978,0.852649424,1.052978981,
                1.307737124,1.692360258,1.746656897,1.806294777,1.872645279,
                1.947699658,2.034515287,2.138158725,2.268008324,2.444794184,
                2.73327664,3.008241985,3.610913003,4.422461221,5.207607281],
                [0.000126256,0.001262563,0.012625968,0.126603479,0.255338522,
                0.388593098,0.529353194,0.681774103,0.852321169,1.052484878,
                1.306951587,1.690924198,1.745105615,1.804609691,1.870802164,
                1.94566582,2.032244498,2.135581318,2.265008934,2.44114961,
                2.728394364,3.001953896,3.600715792,4.405080131,5.181081688],
                [0.00012623,0.001262298,0.012623317,0.126576465,0.255281381,
                0.388499082,0.529210686,0.681564078,0.852011889,1.052019427,
                1.306211802,1.68957244,1.743645451,1.803023705,1.86906763,
                1.943751954,2.030107915,2.133156598,2.262187768,2.437722528,
                2.723805586,2.996046605,3.59114677,4.388796245,5.156271539],
                [0.000126205,0.001262047,0.012620813,0.126550957,0.255227427,
                0.388410317,0.529076149,0.681365824,0.851719984,1.051580207,
                1.305513886,1.688297694,1.742268607,1.801528327,1.867432297,
                1.941947751,2.028093987,2.130871391,2.259529406,2.43449404,
                2.719484627,2.990486565,3.582149695,4.373509077,5.133016208],
                [0.000126181,0.00126181,0.012618445,0.126526833,0.255176401,
                0.388326375,0.528948932,0.681178377,0.851444031,1.05116506,
                1.304854381,1.687093597,1.74096815,1.800116011,1.86588792,
                1.940244051,2.026192447,2.128713996,2.257020154,2.431447397,
                2.715408718,2.985244052,3.573674844,4.359130025,5.111174452],
                [0.000126159,0.001261586,0.012616202,0.126503982,0.255128071,
                0.388246872,0.528828453,0.681000879,0.851182757,1.050772062,
                1.304230204,1.685954461,1.739737897,1.798780024,1.864427117,
                1.938632686,2.024394147,2.126674013,2.254647813,2.428567627,
                2.711557598,2.980292638,3.565678071,4.345580667,5.090621587],
                [0.000126137,0.001261373,0.012614074,0.126482307,0.255082228,
                0.388171466,0.528714192,0.680832557,0.85093502,1.050399485,
                1.303638589,1.684875122,1.738572312,1.797514341,1.863043285,
                1.937106347,2.022690901,2.124742063,2.252401478,2.425841405,
                2.707913179,2.975608749,3.558120081,4.332791341,5.071247134],
                [0.000126117,0.001261171,0.012612053,0.126461718,0.255038686,
                0.388099848,0.52860568,0.680672717,0.850699796,1.050045779,
                1.303077053,1.683851014,1.737466471,1.796313557,1.861730496,
                1.93565848,2.02107537,2.122909811,2.250271372,2.423256774,
                2.704459262,2.971171285,3.55096576,4.320699958,5.052952855],
                [0.000126098,0.001260979,0.012610131,0.126442137,0.254997276,
                0.38803174,0.528502492,0.680520735,0.85047616,1.049709544,
                1.302543359,1.682878003,1.736415819,1.795172806,1.860483419,
                1.934283185,2.019540948,2.121169734,2.248248705,2.420802986,
                2.701181298,2.96696131,3.544183642,4.309251003,5.035651111],
                [0.00012608,0.001260796,0.012608301,0.126423491,0.254957845,
                0.387966889,0.528404246,0.680376045,0.850263276,1.049389519,
                1.302035487,1.681952358,1.735416361,1.794087695,1.85929724,
                1.932975177,2.018081679,2.119515049,2.24632555,2.418470354,
                2.69806618,2.962961776,3.537745445,4.298394685,5.019263469],
                [0.000126062,0.001260622,0.012606556,0.126405715,0.254920254,
                0.387905069,0.528310597,0.680238134,0.850060387,1.049084557,
                1.301551608,1.681070704,1.734464444,1.793054252,1.858167606,
                1.931729572,2.016692173,2.117939623,2.244494741,2.416250123,
                2.695102073,2.959157299,3.531625677,4.288086224,5.003719532],
                [0.000126046,0.001260455,0.01260489,0.126388749,0.254884377,
                0.38784607,0.528221227,0.680106537,0.849866805,1.048793622,
                1.30109006,1.680229977,1.733556755,1.792068872,1.857090568,
                1.930542038,2.015367547,2.116437895,2.242749784,2.414134361,
                2.692278259,2.955533954,3.525801306,4.278285232,4.988955935],
                [0.00012603,0.001260296,0.012603299,0.126372539,0.2548501,
                0.387789704,0.528135852,0.679980831,0.849681905,1.048515765,
                1.300649332,1.679427393,1.732690282,1.79112828,1.856062536,
                1.929408607,2.014103359,2.115004812,2.241084781,2.412115869,
                2.689585012,2.952079114,3.520251464,4.26895519,4.974915484],
                [0.000126014,0.001260144,0.012601777,0.126357036,0.254817318,
                0.3877358,0.528054208,0.679860627,0.849505115,1.048250126,
                1.300228048,1.678660414,1.731862281,1.790229489,1.855080242,
                1.928325664,2.012895567,2.113635774,2.239494359,2.410188089,
                2.687013484,2.948781309,3.514957205,4.260062999,4.961546435],
                [0.000126,0.001259998,0.01260032,0.126342194,0.254785936,
                0.387684198,0.527976057,0.679745573,0.849335912,1.047995915,
                1.299824947,1.677926722,1.731070245,1.789369774,1.854140699,
                1.927289911,2.01174048,2.112326583,2.237973619,2.408345042,
                2.68455561,2.945630052,3.509901282,4.251578583,4.948801855],
                [0.000125986,0.001259859,0.012598924,0.126327973,0.254755865,
                0.387634755,0.527901179,0.679635345,0.84917382,1.047752411,
                1.299438879,1.677224197,1.730311881,1.788546685,1.853241176,
                1.926298328,2.010634722,2.111073398,2.23651808,2.406581265,
                2.682204018,2.942615802,3.505067969,4.243474559,4.936639077],
                [0.000125972,0.001259725,0.012597585,0.126314333,0.254727026,
                0.387587338,0.527829373,0.679529645,0.849018399,1.047518951,
                1.299068785,1.676550893,1.729585086,1.787757851,1.85237917,
                1.925348149,2.009575199,2.109872701,2.235123638,2.40489175,
                2.679951964,2.939729816,3.50044289,4.235725936,4.925019228],
                [0.00012596,0.001259596,0.0125963,0.12630124,0.254699343,
                0.387541825,0.527760453,0.6794282,0.848869245,1.047294927,
                1.298713694,1.675905026,1.72888793,1.787001214,1.851552381,
                1.924436836,2.008559072,2.108721264,2.233786526,2.403271907,
                2.677793261,2.936964083,3.496012881,4.228309866,4.913906815],
                [0.000125947,0.001259473,0.012595065,0.126288662,0.254672749,
                0.387498103,0.527694249,0.679330759,0.848725986,1.047079777,
                1.298372713,1.675284951,1.728218636,1.786274842,1.850758696,
                1.92356205,2.007583728,2.107616115,2.23250328,2.401717513,
                2.675722224,2.934311242,3.491765862,4.221205414,4.90326937],
                [0.000125935,0.001259354,0.012593878,0.126276569,0.254647181,
                0.387456069,0.527630603,0.679237087,0.848588281,1.046872985,
                1.298045016,1.674689154,1.727575567,1.785576954,1.849996165,
                1.92272164,2.006646761,2.10655452,2.231270736,2.400224681,
                2.67373362,2.931764522,3.487690733,4.214393364,4.89307713],
                [0.000125924,0.00125924,0.012592735,0.126264933,0.254622581,
                0.387415628,0.52756937,0.679146972,0.848455812,1.046674072,
                1.297729843,1.674116237,1.726957211,1.784905909,1.84926299,
                1.921913618,2.005745949,2.105533953,2.230085899,2.398789825,
                2.671822625,2.929317682,3.483777272,4.207856049,4.883302766],
                [0.000125913,0.00125913,0.012591635,0.126253729,0.254598894,
                0.387376689,0.527510416,0.679060214,0.848328287,1.046482598,
                1.297426488,1.673564907,1.72636217,1.784260186,1.848557506,
                1.921136147,2.004879275,2.104552082,2.228946075,2.397409633,
                2.669984784,2.926964961,3.480016049,4.201577194,4.873921136],
                [0.000125902,0.001259024,0.012590575,0.126242933,0.254576071,
                0.387339171,0.527453615,0.67897663,0.848205433,1.046298152,
                1.2971343,1.673033966,1.72578915,1.783638379,1.847878176,
                1.920387528,2.004044769,2.103606746,2.227848747,2.39608104,
                2.668215976,2.924701034,3.476398358,4.195541782,4.864909073],
                [0.000125892,0.001258922,0.012589553,0.126232524,0.254554065,
                0.387302997,0.527398851,0.678896047,0.848086998,1.046120354,
                1.296852673,1.672522304,1.725236952,1.783039183,1.84722357,
                1.919666183,2.003240704,2.102695945,2.226791582,2.394801206,
                2.666512385,2.922520967,3.472916138,4.189735939,4.856245195],
                [0.000125882,0.001258823,0.012588567,0.126222481,0.254532833,
                0.387268098,0.527346018,0.678818308,0.847972749,1.045948852,
                1.296581044,1.672028889,1.72470446,1.782461387,1.846592364,
                1.91897065,2.002465444,2.10181782,2.225772416,2.393567496,
                2.664870469,2.920420189,3.469561926,4.184146823,4.847909735],
                [0.000125873,0.001258728,0.012587616,0.126212785,0.254512335,
                0.387234405,0.527295014,0.678743264,0.847862467,1.045783317,
                1.29631889,1.671552763,1.724190638,1.781903866,1.845983324,
                1.918299566,2.001717468,2.100970643,2.224789235,2.392377461,
                2.66328694,2.91839445,3.466328793,4.178762532,4.839884392],
                [0.000125864,0.001258636,0.012586696,0.126203418,0.254492534,
                0.387201858,0.527245746,0.678670778,0.847755949,1.045623442,
                1.296065725,1.671093033,1.723694522,1.78136557,1.845395304,
                1.917651665,2.000995361,2.100152808,2.223840168,2.391228822,
                2.661758738,2.916439804,3.463210303,4.173572019,4.8321522],
                [0.000125855,0.001258547,0.012585807,0.126194364,0.254473395,
                0.3871704,0.527198127,0.678600721,0.847653006,1.045468943,
                1.295821093,1.670648865,1.723215211,1.780845521,1.844827279,
                1.917025767,2.000297804,2.099362815,2.222923468,2.390119457,
                2.660283014,2.914552572,3.460200467,4.168565017,4.824697406],
                [0.000125846,0.001258461,0.012584947,0.126185608,0.254454885,
                0.387139977,0.527152076,0.678532973,0.847553461,1.045319553,
                1.295584571,1.670219484,1.722751866,1.780342808,1.844278161,
                1.916420769,1.999623567,2.098599268,2.222037509,2.389047385,
                2.658857111,2.912729327,3.457293707,4.16373197,4.817505364],
                [0.000125838,0.001258378,0.012584116,0.126177135,0.254436973,
                0.387110537,0.527107517,0.678467421,0.847457148,1.045175022,
                1.295355762,1.669804163,1.722303703,1.779856578,1.843747063,
                1.915835642,1.998971498,2.097860859,2.22118077,2.388010758,
                2.657478549,2.910966871,3.454484819,4.159063975,4.81056244],
                [0.00012583,0.001258298,0.01258331,0.126168931,0.254419631,
                0.387082036,0.527064377,0.678403961,0.847363912,1.045035117,
                1.295134294,1.669402222,1.721869987,1.779386034,1.843233111,
                1.915269423,1.998340522,2.097146368,2.220351829,2.387007846,
                2.656145008,2.909262212,3.451768943,4.154552726,4.803855928],
                [0.000125822,0.00125822,0.01258253,0.126160984,0.254402833,
                0.387054427,0.527022591,0.678342494,0.847273609,1.044899619,
                1.29491982,1.669013026,1.721450031,1.778930427,1.842735489,
                1.914721209,1.997729633,2.096454652,2.219549355,2.386037031,
                2.65485432,2.907612556,3.449141537,4.150190463,4.797373967],
                [0.000125814,0.001258144,0.012581774,0.126153282,0.254386552,
                0.387027671,0.526982096,0.678282928,0.847186101,1.044768324,
                1.294712013,1.668635976,1.72104319,1.77848906,1.842253429,
                1.914190156,1.997137887,2.095784638,2.218772099,2.385096797,
                2.653604451,2.906015283,3.44659835,4.145969932,4.791105479],
                [0.000125807,0.001258071,0.012581041,0.126145814,0.254370766,
                0.387001728,0.526942832,0.678225175,0.847101262,1.044641039,
                1.294510568,1.668270515,1.720648911,1.778061273,1.841786214,
                1.91367547,1.996564396,2.095135322,2.218018893,2.384185721,
                2.652393497,2.904467937,3.444135396,4.14188434,4.785040101],
                [0.0001258,0.001258,0.012580329,0.126138569,0.254355453,
                0.386976561,0.526904745,0.678169155,0.84701897,1.044517583,
                1.294315197,1.667916115,1.720266525,1.77764645,1.841333168,
                1.913176406,1.996008331,2.094505758,2.217288636,2.383302468,
                2.651219666,2.902968213,3.441748941,4.137927322,4.779168128],
                [0.000125793,0.001257931,0.012579639,0.126131538,0.25434059,
                0.386952137,0.526867782,0.67811479,0.846939114,1.044397786,
                1.294125629,1.667572281,1.719895546,1.77724401,1.840893656,
                1.912692263,1.995468907,2.093895059,2.216580296,2.382445783,
                2.650081279,2.901513948,3.439435475,4.134092907,4.773480466],
                [0.000125786,0.001257864,0.012578969,0.126124711,0.25432616,
                0.386928423,0.526831894,0.678062008,0.846861585,1.044281487,
                1.293941609,1.667238549,1.719535472,1.776853407,1.84046708,
                1.912222383,1.99494539,2.09330239,2.215892901,2.381614484,
                2.648976754,2.900103106,3.437191701,4.130375487,4.767968579],
                [0.00012578,0.001257798,0.012578317,0.126118079,0.254312143,
                0.386905388,0.526797036,0.678010741,0.846786285,1.044168536,
                1.293762898,1.66691448,1.719185829,1.776474127,1.840052879,
                1.911766144,1.994437086,2.092726964,2.215225536,2.38080746,
                2.647904603,2.89873377,3.435014519,4.126769791,4.76262445],
                [0.000125774,0.001257735,0.012577685,0.126111635,0.254298521,
                0.386883004,0.526763162,0.677960925,0.846713118,1.04405879,
                1.293589269,1.666599659,1.71884617,1.776105684,1.839650521,
                1.91132296,1.993943341,2.092168039,2.214577337,2.380023664,
                2.646863423,2.897404136,3.432901008,4.123270861,4.757440544],
                [0.000125767,0.001257674,0.01257707,0.12610537,0.254285278,
                0.386861243,0.526730233,0.677912499,0.846641995,1.043952114,
                1.293420507,1.666293697,1.718516073,1.77574762,1.839259505,
                1.91089228,1.993463539,2.091624915,2.213947491,2.379262106,
                2.645851891,2.896112503,3.430848416,4.119874026,4.752409768],
                [0.000125761,0.001257614,0.012576471,0.126099276,0.2542724,
                0.38684008,0.52669821,0.677865405,0.846572832,1.043848382,
                1.293256413,1.665996224,1.718195141,1.775399503,1.838879359,
                1.910473581,1.992997097,2.091096932,2.213335228,2.378521854,
                2.64486876,2.894857264,3.428854147,4.116574884,4.747525441],
                [0.000125756,0.001257556,0.012575889,0.126093348,0.254259869,
                0.38681949,0.526667054,0.67781959,0.846505548,1.043747473,
                1.293096793,1.665706893,1.717882997,1.775060923,1.838509635,
                1.91006637,1.992543466,2.090583466,2.212739824,2.377802026,
                2.643912849,2.893636903,3.426915752,4.113369284,4.742781266],
                [0.00012575,0.001257499,0.012575323,0.126087578,0.254247674,
                0.386799451,0.526636732,0.677775001,0.846440068,1.043649274,
                1.292941469,1.665425374,1.717579285,1.774731495,1.838149911,
                1.909670182,1.992102124,2.090083926,2.212160592,2.377101787,
                2.642983043,2.892449987,3.425030916,4.110253304,4.738171302],
                [0.000125744,0.001257444,0.012574771,0.12608196,0.2542358,
                0.386779941,0.526607211,0.677731591,0.84637632,1.043553676,
                1.292790268,1.665151354,1.717283668,1.774410853,1.837799787,
                1.909284574,1.991672579,2.089597752,2.21159688,2.376420351,
                2.642078309,2.891295159,3.423197448,4.10722324,4.733689937],
                [0.000125739,0.00125739,0.012574234,0.126076489,0.254224236,
                0.386760938,0.526578459,0.677689313,0.846314237,1.043460579,
                1.292643029,1.664884538,1.716995826,1.774098648,1.837458882,
                1.908909128,1.991254363,2.089124452,2.211048075,2.375756968,
                2.641197607,2.890171135,3.421413279,4.104275589,4.729331869],
                [0.000125734,0.001257338,0.012573711,0.126071157,0.254212968,
                0.386742424,0.526550447,0.677648124,0.846253754,1.043369884,
                1.292499597,1.664624645,1.716715457,1.773794554,1.837126839,
                1.908543447,1.990847036,2.088663453,2.210513591,2.375110931,
                2.64034001,2.889076701,3.419676447,4.101407036,4.725092082],
                [0.000125729,0.001257287,0.0125732,0.126065961,0.254201986,
                0.38672438,0.526523146,0.677607981,0.846194811,1.0432815,
                1.292359828,1.66437141,1.716442273,1.773498257,1.836803315,
                1.908187201,1.990450177,2.088214314,2.209992877,2.374481569,
                2.639504623,2.888010702,3.417985094,4.09861444,4.72096583],
                [0.000125724,0.001257237,0.012572703,0.126060896,0.25419128,
                0.386706788,0.526496529,0.677568846,0.846137348,1.043195341,
                1.292223583,1.664124579,1.716176002,1.773209461,1.836487987,
                1.907839944,1.990063387,2.087776582,2.209485406,2.373868245,
                2.638690591,2.886972044,3.416337455,4.095894825,4.716948615],
                [0.000125719,0.001257189,0.012572218,0.126055955,0.254180838,
                0.386689632,0.526470572,0.677530682,0.846081311,1.043111322,
                1.29209073,1.663883913,1.715916384,1.772927886,1.836180548,
                1.90750138,1.989686288,2.087349829,2.20899068,2.373270352,
                2.637897108,2.88595969,3.41473186,4.093245369,4.713036175],
                [0.000125714,0.001257141,0.012571745,0.126051135,0.254170651,
                0.386672895,0.52644525,0.677493452,0.846026648,1.043029366,
                1.291961144,1.663649185,1.715663174,1.772653264,1.835880705,
                1.907171188,1.989318521,2.086933648,2.208508225,2.372687317,
                2.637123405,2.884972653,3.413166718,4.090663392,4.709224467],
                [0.000125709,0.001257095,0.012571283,0.126046431,0.25416071,
                0.386656562,0.52642054,0.677457122,0.845973309,1.042949397,
                1.291834705,1.663420175,1.715416137,1.772385341,1.83558818,
                1.90684906,1.988959743,2.086527649,2.208037589,2.372118591,
                2.636368751,2.884009995,3.411640522,4.08814635,4.705509653],
                [0.000125705,0.00125705,0.012570832,0.12604184,0.254151007,
                0.386640619,0.52639642,0.67742166,0.845921246,1.042871344,
                1.291711301,1.66319668,1.715175049,1.772123873,1.835302708,
                1.906534704,1.988609629,2.086131464,2.207578344,2.371563655,
                2.635632452,2.883070824,3.410151834,4.085691826,4.701888087],
                [0.000125701,0.001257006,0.012570392,0.126037357,0.254141532,
                0.386625052,0.526372869,0.677387037,0.845870414,1.042795139,
                1.291590824,1.6629785,1.7149397,1.771868684,1.835024038,
                1.906227843,1.988267868,2.085744741,2.207130082,2.371022013,
                2.634913847,2.882154291,3.408699291,4.083297519,4.698356304],
                [0.000125696,0.001256963,0.012569962,0.126032978,0.254132278,
                0.386609848,0.526349868,0.677353221,0.845820769,1.042720716,
                1.291473171,1.66276545,1.714709885,1.77161945,1.834751929,
                1.905928211,1.987934166,2.085367145,2.206692412,2.370493194,
                2.634212304,2.881259586,3.407281591,4.080961241,4.694911008],
                [0.000125692,0.001256921,0.012569542,0.1260287,0.254123236,
                0.386594994,0.526327396,0.677320186,0.845772271,1.042648015,
                1.291358243,1.66255735,1.714485413,1.771376012,1.834486152,
                1.905635557,1.987608241,2.084998357,2.206264964,2.369976746,
                2.633527223,2.880385941,3.405897495,4.078680908,4.691549064],
                [0.000125688,0.00125688,0.012569131,0.126024519,0.254114401,
                0.386580479,0.526305437,0.677287904,0.84572488,1.042576975,
                1.291245948,1.66235403,1.714266099,1.771138171,1.834226488,
                1.905349639,1.987289823,2.084638073,2.205847383,2.369472241,
                2.632858032,2.879532619,3.404545824,4.076454537,4.688267483],
                [0.000125684,0.00125684,0.01256873,0.126020433,0.254105765,
                0.38656629,0.526283973,0.67725635,0.84567856,1.042507542,
                1.291136195,1.662155326,1.714051767,1.770905736,1.83397273,
                1.905070227,1.986978657,2.084286001,2.205439333,2.36897927,
                2.632204185,2.878698921,3.403225452,4.074280234,4.685063421],
                [0.00012568,0.001256801,0.012568338,0.126016437,0.254097321,
                0.386552417,0.526262987,0.6772255,0.845633273,1.04243966,
                1.291028899,1.661961085,1.713842249,1.770678524,1.833724678,
                1.904797102,1.986674497,2.083941864,2.20504049,2.368497442,
                2.631565159,2.877884176,3.401935303,4.072156195,4.681934165],
                [0.000125676,0.001256762,0.012567954,0.12601253,0.254089063,
                0.38653885,0.526242464,0.67719533,0.845588986,1.042373279,
                1.290923979,1.661771156,1.713637386,1.770456362,1.833482141,
                1.904530055,1.98637711,2.083605396,2.204650547,2.368026382,
                2.630940457,2.877087746,3.400674352,4.070080697,4.678877126],
                [0.000125672,0.001256725,0.012567579,0.126008707,0.254080984,
                0.386525579,0.526222388,0.677165819,0.845545666,1.042308349,
                1.290821356,1.661585397,1.713437085,1.770239083,1.833244938,
                1.904268883,1.986086272,2.083276345,2.204269209,2.367565735,
                2.630329602,2.87630902,3.39944162,4.068052096,4.675889835],
                [0.000125669,0.001256688,0.012567212,0.126004967,0.25407308,
                0.386512593,0.526202745,0.677136944,0.845503282,1.042244823,
                1.290720956,1.661403674,1.713241078,1.770026527,1.833012895,
                1.904013396,1.985801768,2.082954467,2.203896193,2.367115157,
                2.629732138,2.875547414,3.398236169,4.06606882,4.672969932],
                [0.000125665,0.001256652,0.012566852,0.126001306,0.254065344,
                0.386499884,0.526183521,0.677108686,0.845461803,1.042182656,
                1.290622708,1.661225856,1.713049284,1.769818543,1.832785845,
                1.90376341,1.985523395,2.082639531,2.203531231,2.366674352,
                2.629147631,2.874802371,3.397057105,4.064129369,4.670115165],
                [0.000125662,0.001256617,0.0125665,0.125997723,0.254057771,
                0.386487444,0.526164702,0.677081025,0.845421202,1.042121805,
                1.290526543,1.661051818,1.712861569,1.769614985,1.832563628,
                1.903518749,1.985250956,2.082331315,2.203174065,2.36624295,
                2.628575664,2.874073355,3.395903569,4.062232305,4.667323378],
                [0.000125658,0.001256583,0.012566156,0.125994214,0.254050356,
                0.386475263,0.526146277,0.677053941,0.84538145,1.042062229,
                1.290432395,1.660881441,1.712677804,1.769415712,1.832346092,
                1.903279245,1.984984263,2.082029605,2.202824447,2.365820681,
                2.628015836,2.873359855,3.394774743,4.060376254,4.664592513],
                [0.000125655,0.001256549,0.012565818,0.125990778,0.254043095,
                0.386463334,0.526128233,0.677027419,0.845342521,1.042003888,
                1.290340202,1.660714611,1.712497867,1.769220591,1.832133091,
                1.903044737,1.984723136,2.081734197,2.202482141,2.365407256,
                2.627467767,2.872661381,3.393669841,4.058559902,4.661920596],
                [0.000125652,0.001256516,0.012565488,0.125987412,0.254035981,
                0.386451649,0.526110558,0.677001439,0.84530439,1.041946743,
                1.290249904,1.660551218,1.712321638,1.769029493,1.831924485,
                1.902815069,1.984467404,2.081444897,2.202146919,2.365002401,
                2.626931088,2.871977464,3.39258811,4.056781988,4.659305742],
                [0.000125648,0.001256483,0.012565164,0.125984114,0.254029012,
                0.3864402,0.526093241,0.676975986,0.845267032,1.041890759,
                1.290161442,1.660391157,1.712149004,1.768842296,1.831720138,
                1.902590094,1.9842169,2.081161517,2.201818565,2.364605852,
                2.62640545,2.871307652,3.391528829,4.055041306,4.656746141],
                [0.000125645,0.001256452,0.012564847,0.125980882,0.254022182,
                0.38642898,0.526076271,0.676951043,0.845230425,1.041835901,
                1.290074761,1.660234327,1.711979857,1.76865888,1.831519921,
                1.902369669,1.983971466,2.080883877,2.201496868,2.364217356,
                2.625890514,2.870651515,3.390491307,4.053336698,4.654240062],
                [0.000125642,0.001256421,0.012564536,0.125977714,0.254015488,
                0.386417984,0.526059638,0.676926596,0.845194545,1.041782135,
                1.289989809,1.660080631,1.711814091,1.768479133,1.831323712,
                1.902153657,1.98373095,2.080611804,2.20118163,2.363836671,
                2.625385957,2.870008637,3.389474879,4.051667055,4.651785845],
                [0.000125639,0.00125639,0.012564231,0.125974609,0.254008925,
                0.386407203,0.526043331,0.67690263,0.845159372,1.041729428,
                1.289906533,1.659929976,1.711651607,1.768302946,1.831131389,
                1.901941928,1.983495205,2.080345133,2.200872656,2.363463562,
                2.624891469,2.869378621,3.38847891,4.050031312,4.649381896],
                [0.000125636,0.00125636,0.012563932,0.125971563,0.25400249,
                0.386396632,0.526027343,0.67687913,0.845124884,1.04167775,
                1.289824884,1.659782274,1.711492307,1.768130214,1.83094284,
                1.901734354,1.98326409,2.080083704,2.200569763,2.363097807,
                2.62440675,2.868761084,3.387502788,4.048428448,4.647026688],
                [0.000125633,0.001256331,0.012563639,0.125968577,0.253996179,
                0.386386265,0.526011662,0.676856084,0.845091062,1.04162707,
                1.289744816,1.659637437,1.7113361,1.767960836,1.830757954,
                1.901530816,1.983037471,2.079827365,2.200272772,2.36273919,
                2.623931515,2.86815566,3.386545926,4.046857481,4.644718754],
                [0.00012563,0.001256302,0.012563351,0.125965647,0.253989988,
                0.386376095,0.525996281,0.676833479,0.845057887,1.041577361,
                1.289666283,1.659495384,1.711182895,1.767794715,1.830576625,
                1.901331196,1.982815217,2.079575967,2.199981512,2.362387504,
                2.623465488,2.867561995,3.385607759,4.045317468,4.642456686],
                [0.000125627,0.001256274,0.012563069,0.125962773,0.253983914,
                0.386366118,0.525981191,0.676811301,0.845025341,1.041528594,
                1.289589241,1.659356034,1.711032608,1.767631759,1.830398752,
                1.901135383,1.982597204,2.07932937,2.199695821,2.36204255,
                2.623008403,2.866979751,3.384687744,4.043807505,4.640239129],
                [0.000125625,0.001256246,0.012562792,0.125959952,0.253977954,
                0.386356328,0.525966383,0.676789539,0.844993405,1.041480743,
                1.289513648,1.659219312,1.710885155,1.767471877,1.830224237,
                1.900943268,1.982383312,2.079087437,2.199415538,2.361704136,
                2.622560006,2.866408601,3.383785361,4.042326721,4.638064784],
                [0.000125622,0.001256219,0.01256252,0.125957184,0.253972105,
                0.386346719,0.525951851,0.676768181,0.844962062,1.041433783,
                1.289439464,1.659085144,1.710740458,1.767314985,1.830052985,
                1.900754748,1.982173424,2.078850038,2.199140512,2.361372079,
                2.622120052,2.865848231,3.382900107,4.040874277,4.635932401],
                [0.000125619,0.001256192,0.012562253,0.125954467,0.253966362,
                0.386337287,0.525937585,0.676747216,0.844931297,1.041387688,
                1.289366649,1.658953459,1.710598439,1.767160998,1.829884907,
                1.900569722,1.98196743,2.078617045,2.198870597,2.3610462,
                2.621688304,2.865298338,3.382031498,4.039449368,4.633840776],
                [0.000125617,0.001256166,0.012561991,0.125951799,0.253960725,
                0.386328026,0.52592358,0.676726633,0.844901093,1.041342435,
                1.289295166,1.658824188,1.710459025,1.767009836,1.829719914,
                1.900388095,1.981765221,2.078388338,2.198605651,2.36072633,
                2.621264535,2.864758632,3.381179071,4.038051219,4.631788753],
                [0.000125614,0.00125614,0.012561734,0.125949179,0.253955189,
                0.386318933,0.525909827,0.676706422,0.844871436,1.041298002,
                1.289224979,1.658697266,1.710322146,1.766861423,1.829557922,
                1.900209774,1.981566695,2.078163799,2.198345539,2.360412303,
                2.620848525,2.864228832,3.380342376,4.036679084,4.629775218],
                [0.000125611,0.001256115,0.012561481,0.125946606,0.253949752,
                0.386310003,0.525896321,0.676686574,0.84484231,1.041254366,
                1.289156053,1.658572629,1.710187731,1.766715685,1.829398851,
                1.900034669,1.981371752,2.077943316,2.198090129,2.36010396,
                2.620440064,2.863708667,3.379520983,4.035332242,4.627799099],
                [0.000125609,0.00125609,0.012561233,0.125944079,0.253944411,
                0.38630123,0.525883054,0.676667077,0.844813702,1.041211505,
                1.289088355,1.658450217,1.710055717,1.766572549,1.829242622,
                1.899862694,1.981180296,2.077726779,2.197839295,2.35980115,
                2.620038948,2.863197878,3.378714476,4.034010002,4.625859363],
                [0.000125607,0.001256066,0.012560989,0.125941596,0.253939165,
                0.386292612,0.525870021,0.676647924,0.844785597,1.0411694,
                1.289021851,1.65832997,1.709926038,1.766431946,1.829089159,
                1.899693766,1.980992234,2.077514083,2.197592914,2.359503723,
                2.61964498,2.862696213,3.377922453,4.032711695,4.623955014],
                [0.000125604,0.001256042,0.01256075,0.125939156,0.253934009,
                0.386284144,0.525857214,0.676629104,0.844757982,1.041128031,
                1.28895651,1.658211831,1.709798633,1.766293811,1.82893839,
                1.899527805,1.980807476,2.077305127,2.197350871,2.359211539,
                2.619257971,2.86220343,3.377144529,4.031436677,4.622085093],
                [0.000125602,0.001256018,0.012560514,0.125936759,0.253928943,
                0.386275823,0.525844629,0.67661061,0.844730846,1.041087378,
                1.288892302,1.658095745,1.709673444,1.766158078,1.828790244,
                1.899364732,1.980625937,2.077099814,2.19711305,2.358924459,
                2.618877739,2.861719295,3.376380329,4.030184328,4.620248676],
                [0.0001256,0.001255995,0.012560283,0.125934402,0.253923963,
                0.386267644,0.52583226,0.676592433,0.844704175,1.041047423,
                1.288829199,1.657981659,1.709550412,1.766024685,1.828644654,
                1.899204474,1.980447532,2.076898049,2.196879343,2.358642351,
                2.618504107,2.861243582,3.375629494,4.028954049,4.61844487],
                [0.000125597,0.001255973,0.012560055,0.125932085,0.253919068,
                0.386259603,0.5258201,0.676574565,0.844677957,1.041008148,
                1.288767171,1.657869523,1.709429483,1.765893573,1.828501606,
                1.899046958,1.980272226,2.076699741,2.196649644,2.358365086,
                2.618136904,2.860776074,3.374891676,4.027745261,4.616672815],
                [0.000125595,0.00125595,0.012559832,0.125929808,0.253914256,
                0.386251698,0.525808146,0.676556998,0.844652182,1.040969536,
                1.288706191,1.657759285,1.709310603,1.765764683,1.828360932,
                1.898892115,1.980099853,2.076504801,2.19642385,2.358092542,
                2.617775967,2.86031656,3.374166541,4.026557406,4.61493168],
                [0.000125593,0.001255928,0.012559612,0.125927568,0.253909523,
                0.386243925,0.525796391,0.676539725,0.844626838,1.04093157,
                1.288646234,1.6576509,1.70919372,1.765637959,1.828222625,
                1.898739876,1.979930381,2.076313145,2.196201863,2.357824599,
                2.617421135,2.859864837,3.373453763,4.025389946,4.613220662],
                [0.000125591,0.001255907,0.012559396,0.125925366,0.253904869,
                0.386236281,0.52578483,0.676522738,0.844601914,1.040894235,
                1.288587273,1.65754432,1.709078785,1.765513348,1.828086623,
                1.898590177,1.979763738,2.07612469,2.195983587,2.357561141,
                2.617072256,2.859420708,3.372753029,4.024242359,4.611538988],
                [0.000125588,0.001255885,0.012559183,0.125923199,0.253900291,
                0.386228762,0.52577346,0.67650603,0.8445774,1.040857515,
                1.288529284,1.6574395,1.708965749,1.765390796,1.827952871,
                1.898442955,1.979599854,2.075939356,2.195768929,2.357302056,
                2.616729181,2.858983984,3.372064037,4.023114143,4.609885909],
                [0.000125586,0.001255864,0.012558974,0.125921068,0.253895788,
                0.386221366,0.525762275,0.676489594,0.844553286,1.040821394,
                1.288472243,1.657336398,1.708854565,1.765270254,1.827821312,
                1.898298149,1.97943866,2.075757068,2.195557802,2.357047237,
                2.616391766,2.85855448,3.371386495,4.022004811,4.608260702],
                [0.000125584,0.001255844,0.012558768,0.125918971,0.253891358,
                0.386214089,0.525751271,0.676473425,0.844529562,1.040785858,
                1.288416127,1.657234971,1.708745189,1.765151672,1.827691893,
                1.898155699,1.979280091,2.07557775,2.195350118,2.356796579,
                2.616059873,2.85813202,3.370720118,4.020913893,4.606662669],
                [0.000125582,0.001255824,0.012558565,0.125916908,0.253886998,
                0.386206929,0.525740443,0.676457514,0.844506219,1.040750893,
                1.288360913,1.657135179,1.708637576,1.765035002,1.827564563,
                1.898015549,1.979124084,2.075401331,2.195145793,2.35654998,
                2.615733366,2.857716431,3.370064634,4.019840935,4.605091134],
                [0.00012558,0.001255804,0.012558366,0.125914878,0.253882708,
                0.386199882,0.525729787,0.676441857,0.844483248,1.040716485,
                1.288306581,1.657036982,1.708531685,1.7649202,1.827439271,
                1.897877643,1.978970576,2.07522774,2.194944748,2.356307344,
                2.615412116,2.857307547,3.369419778,4.018785497,4.603545443],
                [0.000125578,0.001255784,0.01255817,0.12591288,0.253878486,
                0.386192947,0.5257193,0.676426447,0.84446064,1.040682621,
                1.288253109,1.656940344,1.708427474,1.764807219,1.827315969,
                1.897741929,1.978819508,2.07505691,2.194746903,2.356068575,
                2.615095998,2.856905207,3.368785292,4.017747153,4.602024964],
                [0.000125576,0.001255765,0.012557976,0.125910912,0.253874329,
                0.38618612,0.525708976,0.676411278,0.844438386,1.040649289,
                1.288200477,1.656845227,1.708324904,1.764696019,1.827194609,
                1.897608353,1.978670823,2.074888776,2.194552183,2.355833582,
                2.614784888,2.856509256,3.36816093,4.016725491,4.600529086],
                [0.000125575,0.001255746,0.012557786,0.125908976,0.253870237,
                0.3861794,0.525698813,0.676396346,0.844416478,1.040616475,
                1.288148665,1.656751594,1.708223936,1.764586555,1.827075146,
                1.897476867,1.978524465,2.074723275,2.194360514,2.355602275,
                2.614478669,2.856119542,3.36754645,4.015720113,4.599057216],
                [0.000125573,0.001255727,0.012557599,0.125907069,0.253866208,
                0.386172782,0.525688807,0.676381643,0.844394908,1.040584169,
                1.288097654,1.656659413,1.708124533,1.764478789,1.826957537,
                1.897347421,1.978380378,2.074560345,2.194171825,2.35537457,
                2.614177227,2.855735919,3.36694162,4.014730632,4.597608781],
                [0.000125571,0.001255709,0.012557415,0.125905191,0.253862241,
                0.386166266,0.525678954,0.676367166,0.844373669,1.040552357,
                1.288047427,1.656568649,1.708026659,1.764372681,1.826841737,
                1.897219968,1.978238512,2.074399926,2.193986048,2.355150381,
                2.61388045,2.855358246,3.366346215,4.013756675,4.596183227],
                [0.000125569,0.00125569,0.012557233,0.125903342,0.253858334,
                0.386159849,0.52566925,0.676352908,0.844352752,1.04052103,
                1.287997964,1.65647927,1.707930279,1.764268193,1.826727706,
                1.897094463,1.978098814,2.074241962,2.193803115,2.354929629,
                2.613588231,2.854986384,3.365760015,4.012797878,4.594780015],
                [0.000125567,0.001255673,0.012557054,0.125901521,0.253854485,
                0.386153529,0.525659692,0.676338865,0.844332151,1.040490175,
                1.287949248,1.656391245,1.707835358,1.764165289,1.826615404,
                1.89697086,1.977961236,2.074086396,2.193622961,2.354712235,
                2.613300466,2.854620202,3.36518281,4.011853892,4.593398624],
                [0.000125565,0.001255655,0.012556878,0.125899727,0.253850694,
                0.386147303,0.525650278,0.676325033,0.844311858,1.040459782,
                1.287901264,1.656304542,1.707741865,1.764063931,1.826504791,
                1.896849119,1.97782573,2.073933174,2.193445524,2.354498123,
                2.613017054,2.854259569,3.364614393,4.010924375,4.59203855],
                [0.000125564,0.001255638,0.012556705,0.125897959,0.25384696,
                0.386141169,0.525641003,0.676311405,0.844291867,1.040429842,
                1.287853994,1.656219133,1.707649767,1.763964087,1.82639583,
                1.896729195,1.977692248,2.073782243,2.193270743,2.354287219,
                2.612737896,2.853904362,3.364054565,4.010008997,4.590699303],
                [0.000125562,0.00125562,0.012556534,0.125896218,0.25384328,
                0.386135126,0.525631865,0.676297979,0.844272171,1.040400344,
                1.287807422,1.656134989,1.707559033,1.763865722,1.826288484,
                1.89661105,1.977560747,2.073633553,2.193098558,2.354079451,
                2.612462898,2.853554457,3.363503134,4.009107438,4.589380408],
                [0.00012556,0.001255604,0.012556365,0.125894502,0.253839654,
                0.386129171,0.52562286,0.676284749,0.844252763,1.040371278,
                1.287761534,1.656052081,1.707469632,1.763768804,1.826182717,
                1.896494644,1.977431183,2.073487053,2.192928913,2.353874751,
                2.612191968,2.853209738,3.362959911,4.008219388,4.588081405],
                [0.000125559,0.001255587,0.012556199,0.125892811,0.253836081,
                0.386123302,0.525613986,0.676271712,0.844233638,1.040342634,
                1.287716314,1.655970383,1.707381537,1.7636733,1.826078495,
                1.896379939,1.977303512,2.073342696,2.19276175,2.35367305,
                2.611925016,2.852870089,3.362424715,4.007344545,4.586801848],
                [0.000125557,0.001255571,0.012556035,0.125891144,0.253832559,
                0.386117518,0.525605241,0.676258862,0.844214788,1.040314405,
                1.287671747,1.655889868,1.707294718,1.763579181,1.825975784,
                1.896266897,1.977177694,2.073200436,2.192597017,2.353474284,
                2.611661954,2.852535401,3.36189737,4.006482618,4.585541303],
                [0.000125555,0.001255555,0.012555874,0.125889501,0.253829088,
                0.386111817,0.52559662,0.676246196,0.844196208,1.04028658,
                1.287627821,1.655810511,1.707209148,1.763486416,1.825874551,
                1.896155484,1.977053689,2.073060226,2.192434696,2.353278389,
                2.611402699,2.852205564,3.361377703,4.005633321,4.584299349],
                [0.000125554,0.001255539,0.012555715,0.125887881,0.253825666,
                0.386106197,0.525588122,0.676233711,0.844177893,1.040259151,
                1.287584521,1.655732288,1.7071248,1.763394976,1.825774765,
                1.896045663,1.976931458,2.072922023,2.192274665,2.353085302,
                2.611147168,2.851880475,3.360865548,4.004796379,4.583075579],
                [0.000125552,0.001255523,0.012555558,0.125886285,0.253822292,
                0.386100656,0.525579743,0.676221401,0.844159836,1.04023211,
                1.287541833,1.655655173,1.707041648,1.763304834,1.825676395,
                1.895937401,1.976810963,2.072785785,2.19211691,2.352894965,
                2.610895282,2.851560031,3.360360744,4.003971526,4.581869596],
                [0.000125551,0.001255507,0.012555404,0.12588471,0.253818965,
                0.386095192,0.525571482,0.676209264,0.844142033,1.040205449,
                1.287499745,1.655579144,1.706959667,1.763215961,1.825579411,
                1.895830665,1.976692167,2.072651468,2.191961383,2.352707318,
                2.610646964,2.851244134,3.359863133,4.0031585,4.580681015],
                [0.000125549,0.001255492,0.012555251,0.125883158,0.253815685,
                0.386089805,0.525563336,0.676197297,0.844124477,1.040179159,
                1.287458244,1.655504178,1.706878832,1.763128331,1.825483785,
                1.895725423,1.976575034,2.072519034,2.191808037,2.352522305,
                2.610402137,2.850932687,3.359372563,4.00235705,4.579509463],
                [0.000125548,0.001255477,0.012555101,0.125881627,0.25381245,
                0.386084492,0.525555303,0.676185494,0.844107165,1.040153233,
                1.287417319,1.655430252,1.70679912,1.763041919,1.825389486,
                1.895621644,1.976459531,2.072388442,2.191656826,2.352339872,
                2.610160729,2.850625598,3.358888886,4.001566929,4.578354576],
                [0.000125546,0.001255462,0.012554953,0.125880116,0.253809259,
                0.386079252,0.52554738,0.676173854,0.84409009,1.040127664,
                1.287376957,1.655357345,1.706720507,1.762956698,1.82529649,
                1.895519298,1.976345623,2.072259655,2.191507706,2.352159963,
                2.60992267,2.850322775,3.358411957,4.0007879,4.577216001],
                [0.000125545,0.001255448,0.012554806,0.125878627,0.253806112,
                0.386074083,0.525539564,0.676162372,0.844073248,1.040102444,
                1.287337146,1.655285437,1.706642971,1.762872645,1.825204767,
                1.895418355,1.976233277,2.072132635,2.191360634,2.351982528,
                2.609687888,2.850024129,3.357941636,4.000019731,4.576093395],
                [0.000125543,0.001255433,0.012554662,0.125877157,0.253803007,
                0.386068984,0.525531855,0.676151046,0.844056635,1.040077566,
                1.287297876,1.655214507,1.706566489,1.762789735,1.825114293,
                1.895318786,1.976122461,2.072007347,2.191215568,2.351807516,
                2.609456318,2.849729577,3.357477786,3.999262195,4.574986424],
                [0.000125542,0.001255419,0.01255452,0.125875708,0.253799944,
                0.386063954,0.525524249,0.676139872,0.844040245,1.040053022,
                1.287259135,1.655144534,1.706491041,1.762707946,1.825025042,
                1.895220564,1.976013145,2.071883755,2.191072468,2.351634877,
                2.609227894,2.849439032,3.357020275,3.998515074,4.573894765],
                [0.000125541,0.001255405,0.012554379,0.125874277,0.253796922,
                0.386058991,0.525516745,0.676128848,0.844024074,1.040028808,
                1.287220914,1.655075501,1.706416605,1.762627255,1.824936989,
                1.895123661,1.975905298,2.071761824,2.190931293,2.351464564,
                2.609002553,2.849152416,3.356568974,3.997778155,4.572818101],
                [0.000125539,0.001255391,0.012554241,0.125872866,0.25379394,
                0.386054094,0.52550934,0.67611797,0.844008118,1.040004915,
                1.2871832,1.655007387,1.706343161,1.76254764,1.824850111,
                1.895028051,1.97579889,2.071641523,2.190792005,2.351296529,
                2.608780231,2.848869647,3.356123758,3.99705123,4.571756127],
                [0.000125538,0.001255378,0.012554104,0.125871473,0.253790997,
                0.386049261,0.525502033,0.676107235,0.843992372,1.039981337,
                1.287145985,1.654940175,1.70627069,1.762469079,1.824764384,
                1.894933709,1.975693894,2.071522817,2.190654566,2.351130728,
                2.60856087,2.84859065,3.355684503,3.996334096,4.570708542],
                [0.000125536,0.001255364,0.012553969,0.125870099,0.253788093,
                0.386044491,0.525494822,0.676096641,0.843976833,1.039958069,
                1.287109259,1.654873847,1.706199173,1.762391552,1.824679786,
                1.894840609,1.975590281,2.071405676,2.19051894,2.350967115,
                2.60834441,2.848315349,3.355251092,3.995626559,4.569675058],
                [0.000125535,0.001255351,0.012553836,0.125868742,0.253785227,
                0.386039784,0.525487704,0.676086184,0.843961496,1.039935104,
                1.287073012,1.654808386,1.70612859,1.762315039,1.824596294,
                1.894748727,1.975488024,2.071290069,2.190385092,2.350805649,
                2.608130794,2.848043672,3.354823408,3.994928425,4.568655391],
                [0.000125534,0.001255338,0.012553704,0.125867403,0.253782397,
                0.386035137,0.525480679,0.676075863,0.843946358,1.039912436,
                1.287037235,1.654743774,1.706058923,1.762239519,1.824513886,
                1.894658039,1.975387096,2.071175966,2.190252986,2.350646287,
                2.607919966,2.847775546,3.354401339,3.99423951,4.567649266],
                [0.000125532,0.001255325,0.012553574,0.125866081,0.253779604,
                0.38603055,0.525473743,0.676065675,0.843931414,1.03989006,
                1.287001918,1.654679996,1.705990155,1.762164974,1.824432543,
                1.894568522,1.975287473,2.071063338,2.190122588,2.350488989,
                2.607711873,2.847510904,3.353984774,3.993559631,4.566656415],
                [0.000125531,0.001255312,0.012553446,0.125864776,0.253776846,
                0.386026022,0.525466897,0.676055617,0.843916661,1.03986797,
                1.286967053,1.654617035,1.705922269,1.762091385,1.824352242,
                1.894480153,1.975189128,2.070952156,2.189993867,2.350333713,
                2.60750646,2.847249678,3.353573607,3.992888613,4.565676578],
                [0.00012553,0.001255299,0.01255332,0.125863487,0.253774124,
                0.386021551,0.525460136,0.676045686,0.843902095,1.03984616,
                1.286932631,1.654554876,1.705855247,1.762018733,1.824272965,
                1.894392911,1.975092037,2.070842393,2.189866789,2.350180423,
                2.607303678,2.846991802,3.353167734,3.992226282,4.564709499],
                [0.000125529,0.001255287,0.012553195,0.125862215,0.253771436,
                0.386017136,0.525453462,0.67603588,0.843887713,1.039824626,
                1.286898644,1.654493503,1.705789074,1.761947,1.824194692,
                1.894306774,1.974996177,2.070734022,2.189741324,2.350029079,
                2.607103475,2.846737211,3.352767052,3.991572472,4.563754932],
                [0.000125527,0.001255274,0.012553071,0.125860958,0.253768781,
                0.386012776,0.52544687,0.676026197,0.843873511,1.039803362,
                1.286865084,1.654432902,1.705723732,1.761876171,1.824117403,
                1.894221722,1.974901524,2.070627016,2.18961744,2.349879645,
                2.606905803,2.846485844,3.352371464,3.990927019,4.562812635],
                [0.000125526,0.001255262,0.012552949,0.125859718,0.253766159,
                0.386008471,0.525440361,0.676016635,0.843859486,1.039782362,
                1.286831942,1.654373058,1.705659207,1.761806226,1.824041081,
                1.894137733,1.974808055,2.07052135,2.189495108,2.349732085,
                2.606710614,2.84623764,3.351980872,3.990289763,4.561882374],
                [0.000125525,0.00125525,0.012552829,0.125858492,0.25376357,
                0.386004219,0.525433932,0.676007191,0.843845635,1.039761623,
                1.286799211,1.654313957,1.705595484,1.761737151,1.823965708,
                1.894054789,1.974715749,2.070416999,2.1893743,2.349586364,
                2.606517862,2.84599254,3.351595182,3.989660551,4.560963919],
                [0.000125524,0.001255238,0.01255271,0.125857282,0.253761012,
                0.386000019,0.525427582,0.675997863,0.843831954,1.039741139,
                1.286766884,1.654255585,1.705532547,1.761668928,1.823891265,
                1.893972869,1.974624584,2.070313938,2.189254987,2.349442448,
                2.606327501,2.845750484,3.351214304,3.98903923,4.560057047],
                [0.000125523,0.001255226,0.012552593,0.125856086,0.253758486,
                0.38599587,0.52542131,0.675988649,0.84381844,1.039720906,
                1.286734952,1.65419793,1.705470382,1.761601542,1.823817736,
                1.893891955,1.974534539,2.070212144,2.189137141,2.349300304,
                2.606139487,2.845511418,3.350838146,3.988425652,4.55916154],
                [0.000125521,0.001255215,0.012552477,0.125854905,0.253755991,
                0.385991772,0.525415114,0.675979547,0.843805091,1.039700918,
                1.286703409,1.654140977,1.705408975,1.761534979,1.823745105,
                1.893812029,1.974445593,2.070111593,2.189020735,2.349159899,
                2.605953776,2.845275286,3.350466622,3.987819675,4.558277186],
                [0.00012552,0.001255203,0.012552362,0.125853738,0.253753525,
                0.385987723,0.525408992,0.675970555,0.843791902,1.039681173,
                1.286672248,1.654084714,1.705348313,1.761469222,1.823673354,
                1.893733072,1.974357726,2.070012263,2.188905744,2.349021201,
                2.605770328,2.845042035,3.350099648,3.987221157,4.557403777],
                [0.000125519,0.001255192,0.012552249,0.125852585,0.253751089,
                0.385983723,0.525402945,0.675961671,0.843778872,1.039661664,
                1.286641461,1.654029129,1.705288381,1.761404257,1.823602467,
                1.893655067,1.974270919,2.069914132,2.188792141,2.348884179,
                2.6055891,2.844811612,3.349737139,3.986629963,4.556541113],
                [0.000125518,0.001255181,0.012552137,0.125851446,0.253748682,
                0.38597977,0.525396969,0.675952892,0.843765998,1.039642389,
                1.286611042,1.653974209,1.705229167,1.761340071,1.82353243,
                1.893577996,1.974185153,2.069817178,2.188679901,2.348748803,
                2.605410053,2.844583965,3.349379015,3.986045957,4.555688996],
                [0.000125517,0.00125517,0.012552027,0.12585032,0.253746304,
                0.385975865,0.525391064,0.675944218,0.843753276,1.039623342,
                1.286580985,1.653919942,1.705170657,1.761276648,1.823463227,
                1.893501844,1.974100409,2.069721379,2.188569001,2.348615044,
                2.605233147,2.844359046,3.349025197,3.985469011,4.554847234],
                [0.000125516,0.001255159,0.012551917,0.125849207,0.253743953,
                0.385972005,0.525385228,0.675935646,0.843740704,1.03960452,
                1.286551283,1.653866318,1.70511284,1.761213977,1.823394843,
                1.893426594,1.974016669,2.069626717,2.188459416,2.348482873,
                2.605058344,2.844136804,3.348675607,3.984898996,4.554015639],
                [0.000125515,0.001255148,0.012551809,0.125848108,0.253741631,
                0.38596819,0.525379461,0.675927175,0.84372828,1.039585919,
                1.286521929,1.653813324,1.705055703,1.761152042,1.823327264,
                1.89335223,1.973933915,2.069533169,2.188351123,2.348352262,
                2.604885608,2.843917194,3.348330171,3.984335788,4.553194028],
                [0.000125514,0.001255137,0.012551703,0.125847021,0.253739335,
                0.38596442,0.525373761,0.675918802,0.843716001,1.039567534,
                1.286492918,1.65376095,1.704999234,1.761090832,1.823260476,
                1.893278736,1.97385213,2.069440718,2.1882441,2.348223183,
                2.604714901,2.843700168,3.347988815,3.983779267,4.552382222],
                [0.000125513,0.001255127,0.012551597,0.125845947,0.253737065,
                0.385960694,0.525368127,0.675910526,0.843703863,1.039549363,
                1.286464244,1.653709184,1.704943421,1.761030334,1.823194465,
                1.893206097,1.973771297,2.069349343,2.188138324,2.34809561,
                2.604546188,2.843485681,3.347651466,3.983229314,4.551580048],
                [0.000125512,0.001255116,0.012551493,0.125844886,0.253734822,
                0.38595701,0.525362558,0.675902346,0.843691866,1.039531402,
                1.286435901,1.653658017,1.704888254,1.760970535,1.823129216,
                1.893134298,1.9736914,2.069259026,2.188033773,2.347969515,
                2.604379435,2.843273688,3.347318055,3.982685814,4.550787334],
                [0.000125511,0.001255106,0.01255139,0.125843836,0.253732605,
                0.385953368,0.525357053,0.675894259,0.843680006,1.039513646,
                1.286407882,1.653607438,1.70483372,1.760911424,1.823064719,
                1.893063326,1.973612422,2.069169749,2.187930426,2.347844875,
                2.604214607,2.843064147,3.346988513,3.982148653,4.550003914],
                [0.00012551,0.001255096,0.012551288,0.125842798,0.253730412,
                0.385949768,0.52535161,0.675886264,0.843668281,1.039496093,
                1.286380184,1.653557436,1.70477981,1.760852988,1.823000958,
                1.892993165,1.973534347,2.069081494,2.187828264,2.347721663,
                2.604051671,2.842857016,3.346662773,3.981617721,4.549229626],
                [0.000125509,0.001255086,0.012551187,0.125841772,0.253728245,
                0.385946209,0.525346229,0.67587836,0.843656689,1.039478738,
                1.286352799,1.653508002,1.704726512,1.760795217,1.822937923,
                1.892923801,1.973457161,2.068994243,2.187727264,2.347599856,
                2.603890595,2.842652252,3.34634077,3.981092911,4.548464311],
                [0.000125508,0.001255076,0.012551088,0.125840758,0.253726102,
                0.38594269,0.525340908,0.675870544,0.843645228,1.03946158,
                1.286325724,1.653459127,1.704673816,1.760738098,1.822875599,
                1.892855223,1.973380848,2.06890798,2.187627408,2.347479429,
                2.603731347,2.842449815,3.34602244,3.980574117,4.547707814],
                [0.000125507,0.001255066,0.012550989,0.125839755,0.253723983,
                0.38593921,0.525335647,0.675862816,0.843633895,1.039444613,
                1.286298952,1.653410801,1.704621713,1.760681621,1.822813977,
                1.892787415,1.973305393,2.068822687,2.187528677,2.347360359,
                2.603573897,2.842249666,3.34570772,3.980061237,4.546959984],
                [0.000125506,0.001255056,0.012550892,0.125838763,0.253721887,
                0.385935769,0.525330445,0.675855175,0.843622688,1.039427836,
                1.286272479,1.653363014,1.704570191,1.760625775,1.822753043,
                1.892720365,1.973230782,2.068738349,2.187431051,2.347242624,
                2.603418213,2.842051767,3.345396549,3.979554169,4.546220672],
                [0.000125505,0.001255047,0.012550796,0.125837782,0.253719815,
                0.385932365,0.5253253,0.675847617,0.843611605,1.039411244,
                1.286246299,1.653315758,1.704519242,1.760570549,1.822692786,
                1.892654061,1.973157001,2.068654949,2.187334512,2.3471262,
                2.603264267,2.84185608,3.345088868,3.979052816,4.545489734],
                [0.000125504,0.001255037,0.0125507,0.125836812,0.253717765,
                0.385929,0.525320212,0.675840143,0.843600644,1.039394835,
                1.286220408,1.653269024,1.704468855,1.760515934,1.822633196,
                1.892588489,1.973084036,2.068572472,2.187239042,2.347011067,
                2.603112029,2.841662567,3.344784619,3.978557082,4.544767028],
                [0.000125503,0.001255028,0.012550606,0.125835853,0.253715738,
                0.385925671,0.525315179,0.675832751,0.843589803,1.039378606,
                1.286194801,1.653222804,1.704419022,1.760461919,1.82257426,
                1.892523639,1.973011873,2.068490903,2.187144623,2.346897203,
                2.602961471,2.841471194,3.344483743,3.978066872,4.544052417],
                [0.000125502,0.001255018,0.012550513,0.125834903,0.253713733,
                0.385922378,0.525310201,0.675825439,0.843579081,1.039362554,
                1.286169474,1.653177088,1.704369734,1.760408494,1.822515969,
                1.892459498,1.9729405,2.068410226,2.187051238,2.346784586,
                2.602812566,2.841281924,3.344186186,3.977582096,4.543345764],
                [0.000125501,0.001255009,0.012550421,0.125833965,0.253711749,
                0.385919121,0.525305277,0.675818207,0.843568474,1.039346676,
                1.286144421,1.653131869,1.704320982,1.76035565,1.822458311,
                1.892396055,1.972869904,2.068330428,2.18695887,2.346673197,
                2.602665286,2.841094722,3.343891892,3.977102662,4.542646938],
                [0.0001255,0.001255,0.012550329,0.125833036,0.253709787,
                0.385915899,0.525300406,0.675811052,0.843557982,1.03933097,
                1.286119638,1.653087139,1.704272756,1.760303437,1.822401278,
                1.892333298,1.972800071,2.068251493,2.186867503,2.346563016,
                2.602519606,2.840909557,3.343600809,3.976628483,4.541955809],
                [0.000125499,0.001254991,0.012550239,0.125832117,0.253707846,
                0.385912712,0.525295588,0.675803974,0.843547603,1.039315431,
                1.286095122,1.653042889,1.704225049,1.760251727,1.822344857,
                1.892271216,1.97273099,2.068173409,2.186777121,2.346454023,
                2.602375498,2.840726393,3.343312883,3.976159474,4.541272251],
                [0.000125498,0.001254982,0.01255015,0.125831208,0.253705925,
                0.385909558,0.52529082,0.675796971,0.843537334,1.039300059,
                1.286070867,1.652999113,1.704177851,1.76020057,1.82228904,
                1.892209799,1.972662649,2.068096161,2.186687707,2.346346199,
                2.602232939,2.840545199,3.343028065,3.97569555,4.54059614],
                [0.000125497,0.001254973,0.012550062,0.125830309,0.253704025,
                0.385906438,0.525286103,0.675790043,0.843527173,1.039284849,
                1.28604687,1.652955802,1.704131156,1.760149957,1.822233817,
                1.892149035,1.972595036,2.068019736,2.186599246,2.346239525,
                2.602091902,2.840365944,3.342746304,3.975236629,4.539927355],
                [0.000125496,0.001254965,0.012549974,0.125829419,0.253702145,
                0.385903351,0.525281436,0.675783188,0.84351712,1.039269801,
                1.286023127,1.65291295,1.704084955,1.76009988,1.822179179,
                1.892088915,1.972528138,2.06794412,2.186511723,2.346133983,
                2.601952364,2.840188597,3.34246755,3.97478263,4.539265778],
                [0.000125496,0.001254956,0.012549888,0.125828538,0.253700285,
                0.385900296,0.525276818,0.675776404,0.843507173,1.03925491,
                1.285999633,1.652870548,1.704039241,1.76005033,1.822125116,
                1.892029428,1.972461946,2.067869302,2.186425123,2.346029555,
                2.6018143,2.840013126,3.342191757,3.974333476,4.538611293],
                [0.000125495,0.001254947,0.012549802,0.125827667,0.253698444,
                0.385897273,0.525272248,0.675769692,0.843497329,1.039240174,
                1.285976384,1.65282859,1.703994004,1.760001299,1.822071619,
                1.891970564,1.972396447,2.067795269,2.186339431,2.345926224,
                2.601677689,2.839839503,3.341918877,3.973889087,4.537963786],
                [0.000125494,0.001254939,0.012549718,0.125826805,0.253696622,
                0.385894281,0.525267725,0.675763049,0.843487588,1.039225592,
                1.285953377,1.652787069,1.703949239,1.759952778,1.82201868,
                1.891912314,1.972331631,2.067722007,2.186254634,2.345823972,
                2.601542506,2.839667698,3.341648864,3.97344939,4.537323147],
                [0.000125493,0.001254931,0.012549634,0.125825951,0.253694819,
                0.38589132,0.525263249,0.675756475,0.843477947,1.039211161,
                1.285930609,1.652745978,1.703904938,1.759904761,1.82196629,
                1.891854668,1.972267488,2.067649506,2.186170718,2.345722782,
                2.60140873,2.839497684,3.341381674,3.97301431,4.536689267],
                [0.000125492,0.001254922,0.012549551,0.125825106,0.253693034,
                0.38588839,0.525258819,0.675749968,0.843468406,1.039196878,
                1.285908074,1.65270531,1.703861094,1.759857238,1.82191444,
                1.891797616,1.972204006,2.067577753,2.186087668,2.345622638,
                2.601276338,2.839329431,3.341117262,3.972583775,4.536062039],
                [0.000125491,0.001254914,0.012549469,0.12582427,0.253691268,
                0.385885489,0.525254434,0.675743528,0.843458962,1.039182741,
                1.285885771,1.65266506,1.703817699,1.759810203,1.821863121,
                1.891741202,1.972141177,2.067506738,2.186005472,2.345523525,
                2.60114531,2.839162914,3.340855585,3.97215772,4.535441359],
                [0.000125491,0.001254906,0.012549387,0.125823443,0.253689519,
                0.385882618,0.525250094,0.675737153,0.843449614,1.039168748,
                1.285863694,1.65262522,1.703774746,1.759763648,1.821812327,
                1.891685313,1.972078988,2.067436447,2.185924116,2.345425425,
                2.601015625,2.838998105,3.340596601,3.971736064,4.534827125],
                [0.00012549,0.001254898,0.012549307,0.125822624,0.253687789,
                0.385879776,0.525245798,0.675730842,0.84344036,1.039154896,
                1.285841842,1.652585784,1.70373223,1.759717565,1.821762048,
                1.891629991,1.972017432,2.067366872,2.185843587,2.345328324,
                2.600887261,2.838834978,3.340340268,3.971318744,4.534219238],
                [0.000125489,0.00125489,0.012549227,0.125821813,0.253686075,
                0.385876963,0.525241545,0.675724596,0.8434312,1.039141185,
                1.285820209,1.652546747,1.703690143,1.759671948,1.821712278,
                1.891575229,1.971956498,2.067298,2.185763874,2.345232207,
                2.600760199,2.838673508,3.340086547,3.970905695,4.533617599]]

        def __init__(self):
                pass

#95% CI
#input:Effect Size(OR,RR,RD),standard error,type of standard error('normal','ln','log')
#output:[95% CI]
def Stat_CI95 (ES, SE, EStype='normal') :
    if EStype=='nomal':
        return [ES-1.96*SE, ES+1.96*SE]
    elif EStype=='ln':
        return [math.exp (math.log(ES)-1.96*SE), math.exp (math.log(ES)+1.96*SE)]
    elif EStype=='log':
        return [math.pow (10, math.log10(ES)-1.96*SE), math.pow (10, math.log10(ES)+1.96*SE)]
    else:
        return [ES-1.96*SE, ES+1.96*SE]
    
#chisquare Test,NOT used in this case
#Input: [Observer values],[Expected values]
#Output: chisquare, k, p
def Stat_chisquare (Obsrd,Expd) :
    idx=0;chisquare=0
    for idx in range(len(Obsrd)) :
        chisquare+=math.pow(Obsrd[idx]-Expd[idx],2)/Expd[idx]
    return [chisquare, len(Obsrd), Stat_guess_CSP(chisquare, len(Obsrd))]

#guess P value from chisquare and k
def Stat_guess_CSP(cs,k):
    lmtbl = lmtbl_chisquare()
    k = int(k)
    k = 2 if k<2 else k
    k = 201 if k>201 else k
    lmt_p_cs=[]
    seeds=lmtbl.seeds
    lmt_p_cs.extend(lmtbl.data[k-2][:])
    idx=0
    if cs<=min(lmt_p_cs):
        return "0.999"
    if cs>=297 :
        return "0.000"
    if cs>=max(lmt_p_cs) :
        return "0.000"
    for idx in range(len(lmt_p_cs)):
        if cs==lmt_p_cs[idx] :
            return '{:.3f}'.format(seeds[idx])
        if cs>lmt_p_cs[idx] :
            if cs<lmt_p_cs[idx+1]:
                p=seeds[idx+1]-((seeds[idx+1]-seeds[idx])*(lmt_p_cs[idx+1]-cs))/(lmt_p_cs[idx+1]-lmt_p_cs[idx])                
                return '{:.3f}'.format(p)

    return '0.000'

#I^2
def Stat_Isquare(Q,k):
    if Q==0:
        return 0
    return max(100*(Q-k+1)/Q,0)

#z test(t/u test)
def Stat_ztest(ES, SE, k=0):
    z=abs(ES/SE)
    p=Stat_guess_TP(z)
    return [z,p]

#guess P value for t test and u test
def Stat_guess_TP(z):
    #from metabase import lmtbl_ttest as lmtbl


#find p-value for two-tailed test

    lmtbl = lmtbl_ttest()
    seeds=lmtbl.seeds
    lmt_p_t=[]
    lmt_p_t=lmtbl.data
    idx=0
    if z<=min(lmt_p_t):
        return "0.999"
    if z>=max(lmt_p_t) :
        p = scipy.stats.norm.sf(abs(z))*2
        # p=seeds[idx+1]-((seeds[idx+1]-seeds[idx])*(lmt_p_t[idx+1]-z))/(lmt_p_t[idx+1]-lmt_p_t[idx])
        return '{:.64f}'.format(p)
    for idx in range(len(lmt_p_t)):
        if z==lmt_p_t[idx] :
            return '{:.8f}'.format(seeds[idx])
        if z>lmt_p_t[idx] :
            if z<lmt_p_t[idx+1]:
                p = scipy.stats.norm.sf(abs(z))*2
                # p=seeds[idx+1]-((seeds[idx+1]-seeds[idx])*(lmt_p_t[idx+1]-z))/(lmt_p_t[idx+1]-lmt_p_t[idx])
                return '{:.64f}'.format(p)

###########################################
'''
For categorical data:
OR, SE(ln(OR));
RR, SE(ln(RR));
RD, SE(RD)
'''

#e1,n1,e2,n2 --> a,b,c,d 
#experiment group:a=e1;b=n1-e1
#control group:c=e2;d=n2-e2
def CATE_en2abcd (studies,type='bad') :
    idx=0;abcd=[]
    for idx in range(len(studies)):
        e1,n1,e2,n2=studies[idx][0:4]
        a,b= (e1,n1-e1) if type=='bad' else (n1-e1,e1)
        c,d= (e2,n2-e2) if type=='bad' else (n2-e2,e2)
        abcd.append([a,b,c,d,studies[idx][4]])
    return abcd

#check empty cell for a,b,c,d
#0.5 should be added to all cells for that study, except when a=c=0 or b=d=0
def CATE_chkzero (study):    
    a,b,c,d=study[0:4]
    if a+b==0 or c+d==0:  #n1=0 or n2=0? should be wrong
        Err="error.(data wrong)"
        raise Exception(Err)
    if (a*b*c*d)==0 :
        if not ((a==0 and c==0) or (b==0 and d==0)):
            a+=0.5;b+=0.5;c+=0.5;d+=0.5
    return [a,b,c,d]

#Odds Ratio
#Input structure:[a,b,c,d]
def CATE_OR (study):
    a,b,c,d=CATE_chkzero (study[0:4])
    return (a*d)/(b*c)

#LnSE of OR
#Input structure:[a,b,c,d]
def CATE_LnSE_OR (study):
    a,b,c,d=CATE_chkzero (study[0:4])
    return math.sqrt(1/a+1/b+1/c+1/d)

#Risk Ratio
#Input structure:[a,b,c,d]
def CATE_RR (study):
    a,b,c,d=CATE_chkzero (study[0:4])
    return (a/(a+b))/(c/(c+d))

#LnSE of Risk Ratio
#Input structure:[a,b,c,d]
def CATE_LnSE_RR (study):
    a,b,c,d=CATE_chkzero (study[0:4])
    return math.sqrt(1/a+1/c-(1/(a+b))-(1/(c+d)))

#Risk Difference
#Input structure:[a,b,c,d]
def CATE_RD (study):
    a,b,c,d=CATE_chkzero (study[0:4])
    return (a/(a+b))-(c/(c+d))

#SE of RD
#Input structure:[a,b,c,d]
def CATE_SE_RD (study):
    a,b,c,d=CATE_chkzero (study[0:4])
    return math.sqrt((a*b/math.pow((a+b),3))+(c*d/math.pow((c+d),3)))

###########################################
'''
For continous data：
SMD, SE(SMD);
MD, SE(MD)
'''
#check empty cell for m1,sd1,n1,m2,sd2,n2
def CONT_chkzero (study):
    #a,b,c,d,e,f=study[0:6]
    if study[2]==0 or study[5]==0:  #n1=0 or n2=0? should be wrong
        return None
    return study

#MD,SE{MD}
#Input structure:[m1,sd1,n1,m2,sd2,n2]
#Output MD,SE{MD}
def CONT_MD (study):
    m1,sd1,n1,m2,sd2,n2=CONT_chkzero (study[0:6])
    N=n1+n2
    md=m1-m2
    se=math.sqrt(sd1*sd1/n1+sd2*sd2/n2)
    return [md,se]

#SMD,SE{SMD}
#Hedges' adjusted g, implemented in RevMan
#Input structure:[m1,sd1,n1,m2,sd2,n2]
#Output SMD,SE{SMD}
def CONT_Heg_SMD (study):
    m1,sd1,n1,m2,sd2,n2=CONT_chkzero (study[0:6])
    # print(m1,sd1,n1,m2,sd2,n2)
    N=n1+n2
    # df=N-1
    df=N-2
    s=math.sqrt(((n1-1)*sd1*sd1+(n2-1)*sd2*sd2)/(N-2))
    md=m1-m2
    # smd=md*(1-3/(4*N-9))/s
    J=(math.gamma(df/2)/(math.sqrt(df/2)*math.gamma((df-1)/2)))
    # smd=md*(math.gamma(df/2)/(math.sqrt(df/2)*math.gamma((df-1)/2)))/s # SMD with gamma function

    # if (s == 0):
    #
    #     #s = 0.0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001
    #     smd = 0
    #     se=(J*J)*math.sqrt( N/(n1*n2)+(smd*smd)/(2*(N-2)) )
    #     return [smd,se]
    if s== 0 :
        s= 0.000000001
    smd=md*J/s # SMD with gamma function

    # se=math.sqrt(N/(n1*n2)+(smd*smd)/(2*(N-3.94))) # des to gia na to vglaleis swsta
    se=(J*J)*math.sqrt( N/(n1*n2)+(smd*smd)/(2*(N-2)) )
    return [smd,se]

#SMD,SE{SMD}
#Cohen's d
#Input structure:[m1,sd1,n1,m2,sd2,n2]
#Output SMD,SE{SMD}
def CONT_Cnd_SMD (study):
    m1,sd1,n1,m2,sd2,n2=CONT_chkzero (study[0:6])
    N=n1+n2
    s=math.sqrt(((n1-1)*sd1*sd1+(n2-1)*sd2*sd2)/(N-2))
    md=m1-m2
    smd=md/s
    se=math.sqrt(N/(n1*n2)+(smd*smd)/(2*(N-2)))
    return [smd,se]

#SMD,SE{SMD}
#Glass's d
#Input structure:[m1,sd1,n1,m2,sd2,n2]
#Output SMD,SE{SMD}
def CONT_Gls_SMD (study):
    m1,sd1,n1,m2,sd2,n2=CONT_chkzero (study[0:6])
    N=n1+n2
    md=m1-m2
    smd=md/sd2
    se=math.sqrt(N/(n1*n2)+(smd*smd)/(2*(n2-1)))
    return [smd,se]

###########################################
'''
M-H
'''
#Functions for M_H method
#Q Test
#input:EFFS=[[es,w],...]:LnOR/LnRR/RD/MD, SE of LnOR/LnRR/RD/MD
#output: Q and p (range string, NOT exactly)
def MH_Q(effs,ttlse):
    Q=0;p=0
    for idx in range(len(effs)):
        Q+=IV_Weight(effs[idx][1])*((effs[idx][0]-ttlse)**2)

    p=Stat_guess_CSP(Q,len(effs))
    return [Q,p]

#====== OR ==============
#Weight of OR
#Input structure:[a,b,c,d]
def MH_Weight_OR (study):
    a,b,c,d=CATE_chkzero (study[0:4])
    return (b*c)/(a+b+c+d)

#Total Effect size & Weight & %95 CI
#Input structure:[[a,b,c,d,'study_name'],...]
#Input models='Fixed': MH;'Random':MH+DL 
#Output structure:[[],...]
#results[0][]:total; [0]'OR', [1]OR, [2]weight, [3]LCI, [4]UCI, [5]N, [6]SE{Ln(OR)}, [7]Q, [8]p for Q test, [9]I^2, [10]z, [11]p for z test...
#results[i][]:studies; [0]study name, [1]OR, [2]weight, [3]LCI, [4]UCI, [5]n, [6]se{Ln(OR)}
def MH_total_OR (studies,models='Fixed'):
    if len(studies)<1:
        Err="error.(no study)"
        raise Exception(Err)
    if len(studies)>200:
        Err="error.(studies should less than 200)"
        raise Exception(Err)
    tw,te,te_r=0,0,0;effs=[]
    w_q,tw_q,te_q=0,0,0    
    r1,r2,r3,r4,r5,r6=0,0,0,0,0,0
    #index=0: "OR",ES,weight,LCI,UCI,N,SE,Q,"p_Q",I^2,z,p_z
    results=[["OR",0,0,0,0,0,0,0,'',0,0,'',None]]
    if models=='Random':
        Tau2=DL_Tausqare (studies,'MH,OR')
        #index=0: "OR",ES,weight,LCI,UCI,N,SE,Q,"p_Q",I^2,z,'p_z',Tau^2
        results=[["OR",0,0,0,0,0,0,0,'',0,0,'',Tau2]]
            
    for idx in range(len(studies)):
        #append a study: study name, effect size, weight, LCI, UCI, n
        es=CATE_OR (studies[idx][0:4])
        se=CATE_LnSE_OR (studies[idx][0:4])
        w=MH_Weight_OR (studies[idx][0:4])
        w_q=w
        if models=='Random':
            w=DL_Weight (se, Tau2)            
        results.append([studies[idx][4]])             #[0]study name
        results[idx+1].append(es)                     #[1]es
        results[idx+1].append(w)                      #[2]weight
        results[idx+1].extend(Stat_CI95 (es,se,'ln')) #[3]LCI,[4]UCI
        results[idx+1].append(sum(studies[idx][0:4])) #[5]n
        results[idx+1].append(se)                     #[6]se{Ln(OR)}
        effs.append([math.log(es),se])                #for cal Q,p
        results[0][2]+=w                              #[0][2]total weight
        results[0][5]+=sum(studies[idx][0:4])         #[0][5]total n
        #for cal total ES
        tw+=w
        te+=es*w
        tw_q+=w_q
        te_q+=es*w_q
        te_r+=math.log(es)*w
        #for cal total SE{ln(OR)}
        a,b,c,d=studies[idx][0:4]
        N=a+b+c+d
        r1+=a*d/N
        r2+=b*c/N
        r3+=((a+d)*a*d/(N*N))
        r4+=((a+d)*b*c/(N*N))
        r5+=((b+c)*a*d/(N*N))
        r6+=((b+c)*b*c/(N*N))
    
    ttlES=te/tw
    ttlSE=math.sqrt((r3/(r1*r1)+(r4+r5)/(r1*r2)+r6/(r2*r2))/2)
    ttlES_Q=te_q/tw_q
    if models=='Random':
        ttlES=math.exp(te_r/tw)
        ttlSE=1/math.sqrt(tw)    
    #[0][1]total ES
    results[0][1]=ttlES
    #[0][3,4]total CI
    results[0][3:5]=Stat_CI95 (ttlES,ttlSE,'ln')
    #[0][6]total SE
    results[0][6]=ttlSE
    #[0][7,8]Q, p of Q_test
    results[0][7:9]=(MH_Q(effs,math.log(ttlES_Q)))
    #[0][9]I^2
    results[0][9]=Stat_Isquare(results[0][7],len(studies))
    #[0][10,11] z, p of z_test    
    results[0][10:12]=(Stat_ztest(math.log(ttlES),ttlSE))
    
    return results

#====== RR ==============
#Weight of RR
#Input structure:[a,b,c,d]
def MH_Weight_RR (study):
    a,b,c,d=CATE_chkzero (study[0:4])
    return ((a+b)*c)/(a+b+c+d)

#Total Effect size & Weight & %95 CI ...
#Input structure:[[a,b,c,d,'study_name'],...]
#Input models='Fixed': MH;'Random':MH+DL 
#Output structure:[[],...]
#results[0][]:total; [0]'RR', [1]RR, [2]weight, [3]LCI, [4]UCI, [5]N, [6]SE{Ln(RR)}, [7]Q, [8]p for Q test, [9]I^2, [10]z, [11]p for z test...
#results[i][]:studies; [0]study name, [1]RR, [2]weight, [3]LCI, [4]UCI, [5]n, [6]se{Ln(RR)}
def MH_total_RR (studies,models='Fixed'):    
    if len(studies)<1:
        Err="error.(no study)"
        raise Exception(Err)
    if len(studies)>200:
        Err="error.(studies should less than 200)"
        raise Exception(Err)
    tw,te,te_r=0,0,0;effs=[]
    w_q,tw_q,te_q=0,0,0   
    r1,r2,r3=0,0,0
    #index=0: "RR",ES,weight,LCI,UCI,N,SE,Q,"p_Q",I^2,z,p_z
    results=[["RR",0,0,0,0,0,0,0,'',0,0,'',None]]
    if models=='Random':
        Tau2=DL_Tausqare (studies,'MH,RR')
        #index=0: "RR",ES,weight,LCI,UCI,N,SE,Q,"p_Q",I^2,z,'p_z',Tau^2
        results=[["RR",0,0,0,0,0,0,0,'',0,0,'',Tau2]]
    
    for idx in range(len(studies)):
        #append a study: study name, effect size, weight, LCI, UCI, n
        es=CATE_RR (studies[idx][0:4])
        se=CATE_LnSE_RR (studies[idx][0:4])
        w=MH_Weight_RR (studies[idx][0:4])
        w_q=w
        if models=='Random':
            w=DL_Weight (se, Tau2)            
        results.append([studies[idx][4]])             #[0]study name
        results[idx+1].append(es)                     #[1]es
        results[idx+1].append(w)                      #[2]weight
        #import pdb;pdb.set_trace()
        results[idx+1].extend(Stat_CI95 (es,se,'ln')) #[3]LCI,[4]UCI
        results[idx+1].append(sum(studies[idx][0:4])) #[5]n
        results[idx+1].append(se)                     #[6]se{Ln(RR)}
        effs.append([math.log(es),se])                #for cal Q,p
        results[0][2]+=w                              #[0][2]total weight
        results[0][5]+=sum(studies[idx][0:4])         #[0][5]total n        
        #for cal total ES
        tw+=w
        te+=es*w
        tw_q+=w_q
        te_q+=es*w_q
        te_r+=math.log(es)*w
        #for cal total SE{ln(RR)}
        a,b,c,d=studies[idx][0:4]
        n1=a+b;n2=c+d;N=n1+n2
        r1+=(n1*n2*(a+c)-a*c*N)/(N*N)
        r2+=a*n2/N
        r3+=c*n1/N
   
    ttlES=te/tw
    ttlSE=math.sqrt(r1/(r2*r3))
    ttlES_Q=te_q/tw_q
    if models=='Random':
        ttlES=math.exp(te_r/tw)
        ttlSE=1/math.sqrt(tw)        
    #[0][1]total ES
    results[0][1]=ttlES
    #[0][3,4]total CI
    results[0][3:5]=Stat_CI95 (ttlES,ttlSE,'ln')
    #[0][6]total SE
    results[0][6]=ttlSE
    #[0][7,8]Q, p of Q_test
    results[0][7:9]=(MH_Q(effs,math.log(ttlES_Q)))
    #[0][9]I^2
    results[0][9]=Stat_Isquare(results[0][7],len(studies))
    #[0][10,11] z, p of z_test    
    results[0][10:12]=(Stat_ztest(math.log(ttlES),ttlSE))
    
    return results

#====== RD ==============
#Weight of RD
#Input structure:[a,b,c,d]
def MH_Weight_RD (study):
    a,b,c,d=CATE_chkzero (study[0:4])
    return (a+b)*(c+d)/(a+b+c+d)

#Total Effect size & Weight & %95 CI
#Input structure:[[a,b,c,d,'study_name'],...]
#Input models='Fixed': MH;'Random':MH+DL 
#Output structure:[[],...]
#results[0][]:total; [0]'RD', [1]RD, [2]weight, [3]LCI, [4]UCI, [5]N, [6]SE{RD}, [7]Q, [8]p for Q test, [9]I^2, [10]z, [11]p for z test...
#results[i][]:studies; [0]study name, [1]RD, [2]weight, [3]LCI, [4]UCI, [5]n, [6]se{RD}
def MH_total_RD (studies,models='Fixed'):
    
    if len(studies)<1:
        Err="error.(no study)"
        raise Exception(Err)
    if len(studies)>200:
        Err="error.(studies should less than 200)"
        raise Exception(Err)
    tw,te=0,0;effs=[]
    w_q,tw_q,te_q=0,0,0 
    r1,r2=0,0
    #index=0: "RD",ES,weight,LCI,UCI,N,SE,Q,"p_Q",I^2,z,p_z
    results=[["RD",0,0,0,0,0,0,0,'',0,0,'',None]]
    if models=='Random':
        Tau2=DL_Tausqare (studies,'MH,RD')
        #index=0: "RD",ES,weight,LCI,UCI,N,SE,Q,"p_Q",I^2,z,'p_z',Tau^2
        results=[["RD",0,0,0,0,0,0,0,'',0,0,'',Tau2]]
    
    for idx in range(len(studies)):
        #append a study: study name, effect size, weight, LCI, UCI, n
        es=CATE_RD (studies[idx][0:4])
        se=CATE_SE_RD (studies[idx][0:4])
        w=MH_Weight_RD (studies[idx][0:4])
        w_q=w
        if models=='Random':
            w=DL_Weight (se, Tau2)            
        results.append([studies[idx][4]])             #[0]study name
        results[idx+1].append(es)                     #[1]es
        results[idx+1].append(w)                      #[2]weight
        #import pdb;pdb.set_trace()
        results[idx+1].extend(Stat_CI95 (es,se))      #[3]LCI,[4]UCI
        results[idx+1].append(sum(studies[idx][0:4])) #[5]n
        results[idx+1].append(se)                     #[6]se{RD}
        effs.append([es,se])                          #for cal Q,p
        results[0][2]+=w                              #[0][2]total weight
        results[0][5]+=sum(studies[idx][0:4])         #[0][5]total n        
        #for cal total ES
        tw+=w
        te+=es*w
        tw_q+=w_q
        te_q+=es*w_q
        #for cal total SE
        a,b,c,d=studies[idx][0:4]
        n1=a+b;n2=c+d;N=n1+n2
        r1+=(a*b*n2*n2*n2+c*d*n1*n1*n1)/(n1*n2*N*N)
        r2+=n1*n2/N
    
    ttlES=te/tw
    ttlSE=math.sqrt(r1/(r2*r2))
    ttlES_Q=te_q/tw_q
    if models=='Random':
        ttlES=te/tw
        ttlSE=1/math.sqrt(tw)
        
    #[0][1]total ES
    results[0][1]=ttlES
    #[0][3,4]total CI
    results[0][3:5]=Stat_CI95 (ttlES,ttlSE)
    #[0][6]total SE
    results[0][6]=ttlSE
    #[0][7,8]Q, p of Q_test
    results[0][7:9]=(MH_Q(effs,ttlES_Q))
    #[0][9]I^2
    results[0][9]=Stat_Isquare(results[0][7],len(studies))
    #[0][10,11] z, p of z_test    
    results[0][10:12] = Stat_ztest(ttlES,ttlSE)

    #import pdb;pdb.set_trace()
    return results

#########################################
'''
Peto
'''
#Functions for Peto method
#====== OR ==============
#Ln Peto OR
#Input structure:[a,b,c,d]
def Peto_LnOR (study):
    a,b,c,d=CATE_chkzero (study[0:4])
    n1=a+b;n2=c+d;N=n1+n2
    z=a-n1*(a+c)/N
    v=n1*n2*(a+c)*(b+d)/(N*N*(N-1))
    return z/v

#Ln Standard error of OR
#Input structure:[a,b,c,d]
def Peto_LnSE_OR (study):
    a,b,c,d=CATE_chkzero (study[0:4])
    n1=a+b;n2=c+d;N=n1+n2
    return math.sqrt(N*N*(N-1)/(n1*n2*(a+c)*(b+d)))

#Hypergeometric variances (Weight) of OR
#Input structure:[a,b,c,d]
def Peto_Weight_OR (study):
    a,b,c,d=CATE_chkzero (study[0:4])
    n1=a+b;n2=c+d;N=n1+n2
    v=n1*n2*(a+c)*(b+d)/(N*N*(N-1))
    return v

#Q Test
#input:EFFS=[[effs,w],...]:LnOR,w; ttles=lnOR
#output: Q and p (range string, NOT exactly)
def Peto_Q(effs,ttles):
    Q=0;p=0
    for idx in range(len(effs)):
        Q+=effs[idx][1]*(effs[idx][0]*effs[idx][0]-ttles*ttles)
        
    p=Stat_guess_CSP(Q,len(effs))
    return [Q,p]

#Total Effect size & Weight & %95 CI
#Input structure:[[a,b,c,d,'study_name'],...]
#Output structure:[[],...]
#results[0][]:total; [0]'OR', [1]OR, [2]weight, [3]LCI, [4]UCI, [5]N, [6]SE{Ln(OR)}, [7]Q, [8]p for Q test, [9]I^2, [10]z, [11]p for z test...
#results[i][]:studies; [0]study name, [1]OR, [2]weight, [3]LCI, [4]UCI, [5]n, [6]se{Ln(OR)}
def Peto_total_OR (studies):    
    if len(studies)<1:
        Err="error.(no study)"
        raise Exception(Err)
    if len(studies)>200:
        Err="error.(studies should less than 200)"
        raise Exception(Err)
    tw,te=0,0;effs=[]
    r1,r2,r3,r4,r5,r6=0,0,0,0,0,0
    #index=0: "OR",ES,weight,LCI,UCI,N,SE,Q,"p_Q",I^2,z,p_z
    results=[["OR",0,0,0,0,0,0,0,'',0,0,'',None]]
    
    for idx in range(len(studies)):
        #append a study: study name, effect size, weight, LCI, UCI, n
        es=Peto_LnOR (studies[idx][0:4])
        se=Peto_LnSE_OR (studies[idx][0:4])
        w=Peto_Weight_OR (studies[idx][0:4])
        results.append([studies[idx][4]])             #[0]study name
        results[idx+1].append(math.exp(es))                     #[1]es
        results[idx+1].append(w)                      #[2]weight
        results[idx+1].extend(Stat_CI95 (math.exp(es),se,'ln')) #[3]LCI,[4]UCI
        results[idx+1].append(sum(studies[idx][0:4])) #[5]n
        results[idx+1].append(se)                     #[6]se{Ln(OR)}
        effs.append([es,w])                           #for cal Q,p
        results[0][2]+=w                              #[0][2]total weight
        results[0][5]+=sum(studies[idx][0:4])         #[0][5]total n        
        #for cal total ES,SE
        tw+=w
        te+=es*w
    
    ttlES=math.exp(te/tw)
    ttlSE=1/math.sqrt(tw)
    #[0][1]total ES
    results[0][1]=ttlES
    #[0][3,4]total CI
    results[0][3:5]=Stat_CI95 (ttlES,ttlSE,'ln')
    #[0][6]total SE
    results[0][6]=ttlSE
    #[0][7,8]Q, p of Q_test
    results[0][7:9]=(Peto_Q(effs,te/tw))
    #[0][9]I^2
    results[0][9]=Stat_Isquare(results[0][7],len(studies))
    #[0][10,11] z, p of z_test    
    results[0][10:12]=(Stat_ztest(math.log(ttlES),ttlSE))
    
    return results

#########################################
'''
Inverse variance
'''
#Weight
#input:se
def IV_Weight(se):
    return 1/(se*se)

#Q Test
#input:EFFS=[[es,se],...]:LnOR/LnRR/RD/MD, SE of LnOR/LnRR/RD/MD
#output: Q and p (range string, NOT exactly)
def IV_Q(effs,ttles):
    Q=0;p=0
    for idx in range(len(effs)):
        Q+=IV_Weight(effs[idx][1])*((effs[idx][0]-ttles)**2)

    p=Stat_guess_CSP(Q,len(effs))
    return [Q,p]

#Total Effect size & Weight & %95 CI
#Input structure:[[a,b,c,d,'study_name'],...]
#Input models='Fixed': IV;'Random':IV+DL 
#Output structure:[[],...]
#results[0][]:total; [0]'OR', [1]OR, [2]weight, [3]LCI, [4]UCI, [5]N, [6]SE{RD}, [7]Q, [8]p for Q test, [9]I^2, [10]z, [11]p for z test...
#results[i][]:studies; [0]study name, [1]OR, [2]weight, [3]LCI, [4]UCI, [5]n, [6]se{RD}
def IV_total_OR (studies,models='Fixed'):
    
    if len(studies)<1:
        Err="error.(no study)"
        raise Exception(Err)
    if len(studies)>200:
        Err="error.(studies should less than 200)"
        raise Exception(Err)
    tw,te=0,0;effs=[]
    w_q,tw_q,te_q=0,0,0    
    #index=0: "OR",ES,weight,LCI,UCI,N,SE,Q,"p_Q",I^2,z,'p_z'
    results=[["OR",0,0,0,0,0,0,0,'',0,0,'',None]]
    if models=='Random':
        Tau2=DL_Tausqare (studies,'IV,OR')
        #index=0: "OR",ES,weight,LCI,UCI,N,SE,Q,"p_Q",I^2,z,'p_z',Tau^2
        results=[["OR",0,0,0,0,0,0,0,'',0,0,'',Tau2]]
    
    for idx in range(len(studies)):
        #append a study: study name, effect size, weight, LCI, UCI, n
        es=CATE_OR (studies[idx][0:4])
        se=CATE_LnSE_OR (studies[idx][0:4])
        w=IV_Weight (se)
        w_q=w
        if models=='Random':
            w=DL_Weight (se, Tau2)            
        results.append([studies[idx][4]])             #[0]study name
        results[idx+1].append(es)                     #[1]es
        results[idx+1].append(w)                      #[2]weight
        #import pdb;pdb.set_trace()
        results[idx+1].extend(Stat_CI95 (es,se,'ln'))      #[3]LCI,[4]UCI
        results[idx+1].append(sum(studies[idx][0:4])) #[5]n
        results[idx+1].append(se)                     #[6]lnse{OR}
        effs.append([math.log(es),se])                #for cal Q,p
        results[0][2]+=w                              #[0][2]total weight
        results[0][5]+=sum(studies[idx][0:4])         #[0][5]total n        
        #for cal total ES
        tw+=w
        te+=math.log(es)*w
        tw_q+=w_q
        te_q+=math.log(es)*w_q
    
    ttlES=math.exp(te/tw)
    ttlSE=1/math.sqrt(tw)
    ttlES_Q=math.exp(te_q/tw_q)
    #[0][1]total ES
    results[0][1]=ttlES
    #[0][3,4]total CI
    results[0][3:5]=Stat_CI95 (ttlES,ttlSE,'ln')
    #[0][6]total SE
    results[0][6]=ttlSE
    #[0][7,8]Q, p of Q_test
    results[0][7:9]=(IV_Q(effs,math.log(ttlES_Q)))
    #[0][9]I^2
    results[0][9]=Stat_Isquare(results[0][7],len(studies))
    #[0][10,11] z, p of z_test    
    results[0][10:12] = Stat_ztest(math.log(ttlES),ttlSE)

    #import pdb;pdb.set_trace()
    return results

#Total Effect size & Weight & %95 CI
#Input structure:[[a,b,c,d,'study_name'],...]
#Input models='Fixed': IV;'Random':IV+DL 
#Output structure:[[],...]
#results[0][]:total; [0]'RR', [1]RR, [2]weight, [3]LCI, [4]UCI, [5]N, [6]SE{RD}, [7]Q, [8]p for Q test, [9]I^2, [10]z, [11]p for z test...
#results[i][]:studies; [0]study name, [1]RR, [2]weight, [3]LCI, [4]UCI, [5]n, [6]se{RD}
def IV_total_RR (studies,models='Fixed'):
    
    if len(studies)<1:
        Err="error.(no study)"
        raise Exception(Err)
    if len(studies)>200:
        Err="error.(studies should less than 200)"
        raise Exception(Err)
    tw,te=0,0;effs=[]
    w_q,tw_q,te_q=0,0,0   
    #index=0: "RR",ES,weight,LCI,UCI,N,SE,Q,"p_Q",I^2,z,'p_z'
    results=[["RR",0,0,0,0,0,0,0,'',0,0,'',None]]
    if models=='Random':
        Tau2=DL_Tausqare (studies,'IV,RR')
        #index=0: "RR",ES,weight,LCI,UCI,N,SE,Q,"p_Q",I^2,z,'p_z',Tau^2
        results=[["RR",0,0,0,0,0,0,0,'',0,0,'',Tau2]]
    
    for idx in range(len(studies)):
        #append a study: study name, effect size, weight, LCI, UCI, n
        es=CATE_RR (studies[idx][0:4])
        se=CATE_LnSE_RR (studies[idx][0:4])
        w=IV_Weight (se)
        w_q=w
        if models=='Random':
            w=DL_Weight (se, Tau2)            
        results.append([studies[idx][4]])             #[0]study name
        results[idx+1].append(es)                     #[1]es
        results[idx+1].append(w)                      #[2]weight
        #import pdb;pdb.set_trace()
        results[idx+1].extend(Stat_CI95 (es,se,'ln'))      #[3]LCI,[4]UCI
        results[idx+1].append(sum(studies[idx][0:4])) #[5]n
        results[idx+1].append(se)                     #[6]lnse{OR}
        effs.append([math.log(es),se])                #for cal Q,p
        results[0][2]+=w                              #[0][2]total weight
        results[0][5]+=sum(studies[idx][0:4])         #[0][5]total n        
        #for cal total ES
        tw+=w
        te+=math.log(es)*w
        tw_q+=w_q
        te_q+=math.log(es)*w_q
    
    ttlES=math.exp(te/tw)
    ttlSE=1/math.sqrt(tw)
    ttlES_Q=math.exp(te_q/tw_q)
    #[0][1]total ES
    results[0][1]=ttlES
    #[0][3,4]total CI
    results[0][3:5]=Stat_CI95 (ttlES,ttlSE,'ln')
    #[0][6]total SE
    results[0][6]=ttlSE
    #[0][7,8]Q, p of Q_test
    results[0][7:9]=(IV_Q(effs,math.log(ttlES_Q)))
    #[0][9]I^2
    results[0][9]=Stat_Isquare(results[0][7],len(studies))
    #[0][10,11] z, p of z_test    
    results[0][10:12] = Stat_ztest(math.log(ttlES),ttlSE)

    #import pdb;pdb.set_trace()
    return results

#Total Effect size & Weight & %95 CI
#Input structure:[[a,b,c,d,'study_name'],...]
#Input models='Fixed': IV;'Random':IV+DL 
#Output structure:[[],...]
#results[0][]:total; [0]'RD', [1]RD, [2]weight, [3]LCI, [4]UCI, [5]N, [6]SE{RD}, [7]Q, [8]p for Q test, [9]I^2, [10]z, [11]p for z test...
#results[i][]:studies; [0]study name, [1]RD, [2]weight, [3]LCI, [4]UCI, [5]n, [6]se{RD}
def IV_total_RD (studies,models='Fixed'):
    
    if len(studies)<1:
        Err="error.(no study)"
        raise Exception(Err)
    if len(studies)>200:
        Err="error.(studies should less than 200)"
        raise Exception(Err)
    tw,te=0,0;effs=[]
    w_q,tw_q,te_q=0,0,0    
    #index=0: "RD",ES,weight,LCI,UCI,N,SE,Q,"p_Q",I^2,z,'p_z'
    results=[["RD",0,0,0,0,0,0,0,'',0,0,'',None]]
    if models=='Random':
        Tau2=DL_Tausqare (studies,'IV,RD')
        #index=0: "RD",ES,weight,LCI,UCI,N,SE,Q,"p_Q",I^2,z,'p_z',Tau^2
        results=[["RD",0,0,0,0,0,0,0,'',0,0,'',Tau2]]
    
    for idx in range(len(studies)):
        #append a study: study name, effect size, weight, LCI, UCI, n
        es=CATE_RD (studies[idx][0:4])
        se=CATE_SE_RD (studies[idx][0:4])
        w=IV_Weight (se)
        w_q=w
        if models=='Random':
            w=DL_Weight (se, Tau2)            
        results.append([studies[idx][4]])             #[0]study name
        results[idx+1].append(es)                     #[1]es
        results[idx+1].append(w)                      #[2]weight
        results[idx+1].extend(Stat_CI95 (es,se))      #[3]LCI,[4]UCI
        results[idx+1].append(sum(studies[idx][0:4])) #[5]n
        results[idx+1].append(se)                     #[6]RD
        effs.append([es,se])                          #for cal Q,p
        results[0][2]+=w                              #[0][2]total weight
        results[0][5]+=sum(studies[idx][0:4])         #[0][5]total n        
        #for cal total ES
        tw+=w
        te+=es*w
        tw_q+=w_q
        te_q+=es*w_q
    
    ttlES=te/tw
    ttlSE=1/math.sqrt(tw)
    ttlES_Q=te_q/tw_q
    #[0][1]total ES
    results[0][1]=ttlES
    #[0][3,4]total CI
    results[0][3:5]=Stat_CI95 (ttlES,ttlSE)
    #[0][6]total SE
    results[0][6]=ttlSE
    #[0][7,8]Q, p of Q_test
    results[0][7:9]=(IV_Q(effs,ttlES_Q))
    #[0][9]I^2
    results[0][9]=Stat_Isquare(results[0][7],len(studies))
    #[0][10,11] z, p of z_test    
    results[0][10:12] = Stat_ztest(ttlES,ttlSE)

    #import pdb;pdb.set_trace()
    return results

#Total Effect size & Weight & %95 CI
#Input structure:[[m1,sd1,n1m2,sd2,n2,'study_name'],...]
#Input models='Fixed': IV;'Random':IV+DL 
#Output structure:[[],...]
#results[0][]:total; [0]'MD', [1]MD, [2]weight, [3]LCI, [4]UCI, [5]N, [6]SE{RD}, [7]Q, [8]p for Q test, [9]I^2, [10]z, [11]p for z test...
#results[i][]:studies; [0]study name, [1]MD, [2]weight, [3]LCI, [4]UCI, [5]n, [6]se{RD}
def IV_total_MD (studies,models='Fixed'):
    
    if len(studies)<1:
        Err="error.(no study)"
        raise Exception(Err)
    if len(studies)>200:
        Err="error.(studies should less than 200)"
        raise Exception(Err)
    tw,te=0,0;effs=[]
    w_q,tw_q,te_q=0,0,0    
    #index=0: "MD",ES,weight,LCI,UCI,N,SE,Q,"p_Q",I^2,z,'p_z'
    results=[["MD",0,0,0,0,0,0,0,'',0,0,'',None]]
    if models=='Random':
        Tau2=DL_Tausqare (studies,'IV,MD')
        #index=0: "MD",ES,weight,LCI,UCI,N,SE,Q,"p_Q",I^2,z,'p_z',Tau^2
        results=[["MD",0,0,0,0,0,0,0,'',0,0,'',Tau2]]
    
    for idx in range(len(studies)):
        #append a study: study name, effect size, weight, LCI, UCI, n
        es,se=CONT_MD (studies[idx][0:6])
        w=IV_Weight (se)
        w_q=w
        if models=='Random':
            w=DL_Weight (se, Tau2)            
        results.append([studies[idx][6]])             #[0]study name
        results[idx+1].append(es)                     #[1]es
        results[idx+1].append(w)                      #[2]weight
        results[idx+1].extend(Stat_CI95 (es,se))      #[3]LCI,[4]UCI
        results[idx+1].append(studies[idx][2]+studies[idx][5]) #[5]n
        results[idx+1].append(se)                     #[6]MD
        effs.append([es,se])                          #for cal Q,p
        results[0][2]+=w                              #[0][2]total weight
        results[0][5]+=(studies[idx][2]+studies[idx][5])         #[0][5]total n        
        #for cal total ES
        tw+=w
        te+=es*w
        tw_q+=w_q
        te_q+=es*w_q
    
    ttlES=te/tw
    ttlSE=1/math.sqrt(tw)
    ttlES_Q=te_q/tw_q
    #[0][1]total ES
    results[0][1]=ttlES
    #[0][3,4]total CI
    results[0][3:5]=Stat_CI95 (ttlES,ttlSE)
    #[0][6]total SE
    results[0][6]=ttlSE
    #[0][7,8]Q, p of Q_test
    results[0][7:9]=(IV_Q(effs,ttlES_Q))
    #[0][9]I^2
    results[0][9]=Stat_Isquare(results[0][7],len(studies))
    #[0][10,11] z, p of z_test    
    results[0][10:12] = Stat_ztest(ttlES,ttlSE)

    #import pdb;pdb.set_trace()
    return results

#Total Effect size & Weight & %95 CI
#Input structure:[[m1,sd1,n1m2,sd2,n2,'study_name'],...]
#Input models='Fixed'(Default): IV;'Random':IV+DL
#Input algorithm='Cnd':Cohen's d; 'Heg':Hedges' adjust g (Default); 'Gls':Glass's D
#Output structure:[[],...]
#results[0][]:total; [0]'SMD', [1]SMD, [2]weight, [3]LCI, [4]UCI, [5]N, [6]SE{RD}, [7]Q, [8]p for Q test, [9]I^2, [10]z, [11]p for z test...
#results[i][]:studies; [0]study name, [1]SMD, [2]weight, [3]LCI, [4]UCI, [5]n, [6]se{RD}
def IV_total_SMD (studies,models='Fixed',algo='Heg'):
    
    if len(studies)<1:
        Err="error.(no study)"
        raise Exception(Err)
    if len(studies)>200:
        Err="error.(studies should less than 200)"
        raise Exception(Err)
    tw,te=0,0;effs=[]
    w_q,tw_q,te_q=0,0,0    
    #index=0: "SMD",ES,weight,LCI,UCI,N,SE,Q,"p_Q",I^2,z,'p_z'
    results=[["SMD",0,0,0,0,0,0,0,'',0,0,'',None]]
    if models=='Random':
        Tau2=DL_Tausqare (studies,'IV,SMD'+','+algo)
        #index=0: "SMD",ES,weight,LCI,UCI,N,SE,Q,"p_Q",I^2,z,'p_z',Tau^2
        results=[["SMD",0,0,0,0,0,0,0,'',0,0,'',Tau2]]
    
    for idx in range(len(studies)):
        #append a study: study name, effect size, weight, LCI, UCI, n
        if algo=='Heg':
            es,se=CONT_Heg_SMD (studies[idx][0:6])
        elif algo=='Cnd':
            es,se=CONT_Cnd_SMD (studies[idx][0:6])
        elif algo=='Gls':
            es,se=CONT_Gls_SMD (studies[idx][0:6])
        w=IV_Weight (se)
        w_q=w
        if models=='Random':
            w=DL_Weight (se, Tau2)            
        results.append([studies[idx][6]])             #[0]study name
        results[idx+1].append(es)                     #[1]es
        results[idx+1].append(w)                      #[2]weight
        results[idx+1].extend(Stat_CI95 (es,se))      #[3]LCI,[4]UCI
        results[idx+1].append(studies[idx][2]+studies[idx][5]) #[5]n
        results[idx+1].append(se)                     #[6]MD
        effs.append([es,se])                          #for cal Q,p
        results[0][2]+=w                              #[0][2]total weight
        results[0][5]+=(studies[idx][2]+studies[idx][5])         #[0][5]total n        
        #for cal total ES
        tw+=w
        te+=es*w
        tw_q+=w_q
        te_q+=es*w_q
    
    ttlES=te/tw
    ttlSE=1/math.sqrt(tw)
    ttlES_Q=te_q/tw_q
    #[0][1]total ES
    results[0][1]=ttlES
    #[0][3,4]total CI
    results[0][3:5]=Stat_CI95 (ttlES,ttlSE)
    #[0][6]total SE
    results[0][6]=ttlSE
    #[0][7,8]Q, p of Q_test
    results[0][7:9]=(IV_Q(effs,ttlES_Q))
    #[0][9]I^2
    results[0][9]=Stat_Isquare(results[0][7],len(studies))
    #[0][10,11] z, p of z_test    
    results[0][10:12] = Stat_ztest(ttlES,ttlSE)

    #import pdb;pdb.set_trace()
    return results

#########################################
'''
DerSimonnian Laird
'''
#Functions for D-L method
#Tau^2
#input:studies
def DL_Tausqare(studies,type):
    if type=='MH,OR':
        rults=MH_total_OR(studies)
    elif type=='MH,RR':
        rults=MH_total_RR(studies)
    elif type=='MH,RD':
        rults=MH_total_RD(studies)
    elif type=='IV,OR':
        rults=IV_total_OR(studies)
    elif type=='IV,RR':
        rults=IV_total_RR(studies)
    elif type=='IV,RD':
        rults=IV_total_RD(studies)
    elif type=='IV,MD':
        rults=IV_total_MD (studies)
    elif type=='IV,SMD' or type=='IV,SMD,Heg':
        rults=IV_total_SMD (studies)
    elif type=='IV,SMD,Cnd':
        rults=IV_total_SMD (studies,'Fixed','Cnd')
    elif type=='IV,SMD,Gls':
        rults=IV_total_SMD (studies,'Fixed','Gls')
        
    #rults[i][6]:se
    #rults[0][7]:Q
    Q=0;w=0;ww=0
    for i in range(1,len(rults)):
        tmp_w=1/(rults[i][6]**2)
        w+=tmp_w
        ww+=tmp_w**2

    Q=rults[0][7]
    k=len(studies)
    if Q<=k-1 :
        return 0
    if w-ww/w == 0 :
        return 0
    return max((Q-k+1)/(w-ww/w),0)

#Weight
#input:se,Tau^2
def DL_Weight(se,Tau2):
    return 1/(se*se+Tau2)

#====== OR,RR,RD ==============
#same to MH and IV with 'Random' models


#########################################
'''
Drawing:forest, funnel
'''
#Method of Drawing ForestPlot
#input list structure:es_w_ci
#results[0]:total; [0]'Total', [1]OR, [2]weight, [3]LCI, [4]UCI, [5]N, [6]SE{Ln(OR)}, [7]Q, [8]p, [9]I2 ...
#results[i]:studies; [0]study name, [1]effect size, [2]weight, [3]LCI, [4]UCI, [5]n,[6]quality ...
def Fig_Forest (size,dpi,es_w_ci, titletxt="Meta-analysis Results",no_ttl=False):
    if es_w_ci[0][0] in "OR,RR" :
        def _x_tran0(x):
            return math.log(x)
        def _x_tran1(x):
            return math.exp(x)
    elif es_w_ci[0][0] in "RD,MD,SMD" :
        def _x_tran0(x):
            return x
        def _x_tran1(x):
            return x
    else:
        Err="error.(failed to get effect size while drawing forest plot)"
        raise Exception(Err)

    myfig = plt.figure(linewidth=1, figsize=size, dpi=dpi)  #Frameon=False, num="Forest plot by PythonMeta", 
    myfig.set_size_inches(size)
    plt.title(titletxt)
    ax = gca()
    ax.spines['right'].set_color('none')
    ax.spines['left'].set_color('none')
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    plt.xticks(fontsize=9)
    plt.yticks(fontsize=9)
    xlim=[];y_k=0;subgrp=[]  
    for i in range(len(es_w_ci)):
        xlim.append(es_w_ci[i][3])
        xlim.append(es_w_ci[i][4])
        stdname=es_w_ci[i][0]
        if stdname[0:5]=="<sub>":  #this line is a subgroup
            subgrp.append(es_w_ci[i])

    xmin= _x_tran0(min(xlim))
    xmax= _x_tran0(max(xlim))
    xmax=max(abs(xmin),abs(xmax))
    xmin=-xmax                    
    plt.xlim(xmin*1.1,xmax*1.1)
    ylabel=[i[0].replace("<sub>","") for i in es_w_ci[1:]]
    if no_ttl==True :
        ax.set_yticks(range(len(es_w_ci)))
        ymax=len(es_w_ci)
        ylabel.extend([""])
        y_k=0
    else:
        ax.set_yticks(range(len(es_w_ci)+3))
        ymax=len(es_w_ci)+3
        ylabel.extend(["","Overall","",""])
        y_k=3
    plt.ylim(0, ymax)
    ylabel.reverse()
    ax.set_yticklabels(ylabel)
    ax.set_xticklabels([round(_x_tran1(x),2) for x in ax.get_xticks()])
    plt.plot([0,0], [0,len(es_w_ci)+3], 'black')     

    if len(subgrp)>0:
        weight_all=subgrp[0][2]
    else:
        weight_all=es_w_ci[0][2]
    N=es_w_ci[0][5];k=0;i_subgrp=0
    for i in range(1,len(es_w_ci)):
        stdname=es_w_ci[i][0]
        if stdname[0:5]=="<sub>":  #this line is a subgroup
            i_subgrp+=1
            if i_subgrp>len(subgrp)-1:
                i_subgrp=len(subgrp)-1
            weight_all=subgrp[i_subgrp][2]
            x=[_x_tran0(es_w_ci[i][1]),
               _x_tran0(es_w_ci[i][3]),
               _x_tran0(es_w_ci[i][1]),
               _x_tran0(es_w_ci[i][4]),
               _x_tran0(es_w_ci[i][1])]
            y=[len(es_w_ci)-i+y_k+0.2,
               len(es_w_ci)-i+y_k,
               len(es_w_ci)-i+y_k-0.2,
               len(es_w_ci)-i+y_k,
               len(es_w_ci)-i+y_k+0.2]
            if (es_w_ci[i][9]<50) :
                plt.fill(x,y,color="blue", lw=1)  #filled: I2<50
            else :
                plt.plot(x,y, 'blue', lw=1)       #empty: I2>50

            continue

        #weight
        weight=es_w_ci[i][2]/weight_all

        #shadow X line
        lncolor,lnstyle=("blue","-")
        plt.plot([_x_tran0(es_w_ci[i][3]),_x_tran0(es_w_ci[i][4])], [len(es_w_ci)-i+y_k,len(es_w_ci)-i+y_k], lncolor, linestyle=lnstyle, lw=1)
        #central block
        k=weight*0.2+0.05 
        x=[_x_tran0(es_w_ci[i][1])-k*(xmax*2.2/ymax),
           _x_tran0(es_w_ci[i][1])+k*(xmax*2.2/ymax),
           _x_tran0(es_w_ci[i][1])+k*(xmax*2.2/ymax),
           _x_tran0(es_w_ci[i][1])-k*(xmax*2.2/ymax),
           _x_tran0(es_w_ci[i][1])-k*(xmax*2.2/ymax)]
        y=[len(es_w_ci)-i+y_k+k,
           len(es_w_ci)-i+y_k+k,
           len(es_w_ci)-i+y_k-k,
           len(es_w_ci)-i+y_k-k,
           len(es_w_ci)-i+y_k+k]
        plt.fill(x,y,color=lncolor, lw=1)  #filled: 

    if no_ttl==True:
        pass
    else:
        #draw total ES from es_w_ci[0]
        x=[_x_tran0(es_w_ci[0][1]),
           _x_tran0(es_w_ci[0][3]),
           _x_tran0(es_w_ci[0][1]),
           _x_tran0(es_w_ci[0][4]),
           _x_tran0(es_w_ci[0][1])]
        y=[2.3,2,1.7,2,2.3]
        if (es_w_ci[0][9]<50) :
            plt.fill(x,y,color="black", lw=1)  #filled: I2<50 
        else :
            plt.plot(x,y, 'black', lw=1)       #empty: I2>50 

    plt.xlabel("Favours Experiment  Favours Control       ",fontsize=10) #effect direction,see Cochrane rules
    #plt.ylabel("Studies")
    myfig.tight_layout() 
    return myfig

#Method of Drawing FunnelPlot
#input list structure:[[Effect Size,SE(ln(effs))],...]
#effs[0][0:2]:total
#effs[idx][0:2]:studies
def Fig_Funnel (size,dpi,effs):
    myfig = plt.figure(linewidth=1,figsize=size,dpi=dpi)  #num="Funnel plot",
    myfig.set_size_inches(size)
    x=[];y=[]
    for i in range(1,len(effs)):
        x.append(effs[i][1])
        y.append(effs[i][6])
    lbl,=plt.plot(x,y,"o", lw=1)

    #plt.xlim(min(x),max(x))
    ymax=max(y)+0.2
    plt.ylim(ymax,0)
    plt.plot ([effs[0][1],effs[0][1]],[0,ymax],color="blue", linestyle="--", lw=1)
    plt.plot ([effs[0][1],effs[0][1]-1.96*ymax],[0,ymax],color="blue", linestyle="--", lw=1)
    plt.plot ([effs[0][1],effs[0][1]+1.96*ymax],[0,ymax],color="blue", linestyle="--", lw=1)

    ax = gca()
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')

    plt.xlabel("Effect Size")
    plt.ylabel("Standard Error")
    myfig.tight_layout() 
    return myfig

if __name__ == '__main__':
    pass
