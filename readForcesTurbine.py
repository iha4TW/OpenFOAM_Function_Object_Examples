# Post Processing file for turbine cases.
# Iván Hernández Alayeto.
# MSc student.
# KEIO University & Technical University of Madrid.
# Contact: ivan.hernandez.alayeto@keio.jp OR ivan.hernandez.alayeto@alumnos.upm.es

import pandas as pd
import numpy as np
from pathlib import Path
from scipy import integrate, interpolate

header=['Time', 'Force pressure x','Force pressure y','Force pressure z', 'Force viscous x','Force viscous y','Force viscous z','Force porous x','Force porous y','Force porous z', 'Torque pressure x','Torque pressure y','Torque pressure z', 'Torque viscous x','Torque viscous y','Torque viscous z','Torque porous x','Torque porous y','Torque porous z']

def bindAngle(theta):
    while theta>=360.:
        theta=theta-360.
    while theta<0:
        theta=theta+360.
    return theta
    
class forcesFile:
    #Class containing the forces as a dataframe, the path of the file and the blade to which it belongs
    filePath=Path('Not initialized')
    forces=pd.Dataframe()
    blade=np.NaN
    startTime=np.NaN
    endTime=np.NaN
    def __init__(self,path):
        self.filePath=path
        #Get the character after forcesBlade in the path of the file.
        self.blade=int(self.filePath[self.filePath('forcesBlade')+len('forcesBlade')])
        self.forces=pd.read_csv(self.path, header=None, skiprows=range(3), delimiter="\s+|\(+|\)+",engine='python')
        #As the parsing introduces NaN when reading '(' or ')' scan the columns and remove them.
        toDrop=[]
        for column in self.forces: 
              if self.forces[column].isnull().values.any(): 
                  toDrop.append(column)
        self.forces=self.forces.drop(toDrop,axis=1)
        self.forces.columns=header
        self.startTime=self.forces.Time[0]
        self.endTime=self.forces.Time[len(self.forces)]

class bladeForces:
    #Class that extracts and orders data for a specified blade from the list of all data.
    number=np.NaN
    omega=7.4
    forces=pd.Dataframe(columns=header)
    avgFx=np.NaN
    avgTz=np.NaN
    P=[]
    T=2*np.pi/omega
    def __init__(self,number,forcesFiles):
        self.number=int(number)
        wrong=False
        #Check that arguments are right.
        for i in forcesFiles:
            if i.blade!= self.number:
                wrong=True
        if wrong:
            print('Arguments initializing bladeForces '+str(self.number)+' are wrong. The object will not be initialized.')
        else:
            #Create a dataframe containing things from all runs.
            forcesFiles=pd.DataFrame(forcesFiles,columns=['Objects'])
            #Get start, end and time span, and blade.
            forcesFiles['startTime','endTime','timeSpan','blade']=np.nan
            forcesFiles.startTime=forcesFiles.Objects.apply(lambda x: x.startTime)
            forcesFiles.endTime=forcesFiles.Objects.apply(lambda x: x.endTime)
            forcesFiles.blade=forcesFiles.Objects.apply(lambda x: x.blade)
            forcesFiles.timeSpan=forcesFiles.endTime-forcesFiles.startTime
            #Remove the files not corresponding to this blade.
            forcesFiles=forcesFiles[forcesFiles.blade==self.number]
            #Get different timesteps.
            uniqueStartingTimes=forcesFiles.startTime.unique()
            uniqueStartingTimes.sort()
            rightIdx=[]
            for i in uniqueStartingTimes:
                rightIdx.append(forcesFiles[forcesFiles.startTime==i].idxmax().timeSpan) #Get the index of the entriy with the same starting time and maximum timespan
            #Now, append all data cutting the previous table to the starting of the next.
            for i in range(len(rightIdx-1)):
                rows=forcesFiles.iloc[rightIdx[i]].Objects.forces.Time<forcesFiles.iloc[rightIdx[i+1]].startTime #List of bools                
                self.forces.append(forcesFiles.iloc[rightIdx[i]].Objects.forces[rows], ignore_index=True)
            #And attach the last file.
            self.forces.append(forcesFiles.iloc[rightIdx[len(rightIdx)-1]].Objects.forces)
            
    def calculateAvgFx(self):
        self.avgFx=self.forces['Force pressure x'].mean()+self.forces['Force viscous x'].mean()
        
    def calculateAvgTz(self):
        self.avTzy=self.forces['Torque pressure z'].mean()+self.forces['Force viscous z'].mean()
        
    def calculatePower(self):
        torque=self.forces['Time','Torque pressure z', 'Torque viscous z'].copy()
        torque['Torque']=torque['Torque pressure z']+torque['Torque viscous z']
        torque.drop(['Torque pressure z', 'Torque viscous z'], axis=1)
        #Iterate over all revolutions.
        t=0
        while t<torque.Time[len(torque)-1]:
            rowsThisRevolution=(torque.Time>t) & (torque.Time<(t+self.T))
            self.P.append( integrate.simps( torque.Time[rowsThisRevolution], torque.Torque(rowsThisRevolution) )*self.omega )
            t=t+self.T
        self.P=np.array(self.P)
        
    def calculateAngle(self):
        self.forces['Angle']=self.forces.Time.apply(lambda t: bindAngle( t/self.T*360. ) )

class turbineData:
    #Class that stores the data of each turbine, gathers power and 
    kind=''
    position=''
    def __init__(self, threeBlades):
        #Get the path of the file.
        casePath=Path('.')
        if 'Single' in casePath.absolute():
            self.kind='Single'
        elif 'Forward' in casePath.absolute():
            self.kind='Forward'
        elif 'Backward' in casePath.absolute():
            self.kind='Backward'
        print('Case recognized as '+self.kind)
        
        if self.kind == 'Single':
            self.position='Centered'
            correctAngle=[0,120,240]
            
        elif self.kind == 'Backward':
            if threeBlades[0].number in [1,2,3]:
                self.position='Top'
                correctAngle=[0,120,240]
            else:
                self.position= 'Bottom'
                correctAngle=[0,120,240]
                
        elif self.kind == 'Forward':
            if threeBlades[0].number in [1,2,3]:
                self.position='Top'
                correctAngle=[0,120,240]
            else:
                self.position= 'Bottom'
                correctAngle=[0,120,240]


#Get all the forces files.
thisCase=Path('.')
allForcesFiles=[x for x in thisCase.glob('**/*forces*.dat')]

for file in allForcesFiles:
    allForcesFiles.append(forcesFile(file))

#Search all the blades:
allBlades=np.unique([x.blade for x in allForcesFiles ])
blades=[]

#Create array with dataFrames for each blade.
for i in allBlades:
    blades.append(bladeForces(i,allForcesFiles))
    
blades=np.array(blades)




