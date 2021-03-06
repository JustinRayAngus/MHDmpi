function[varData]=loadData1DVec(dataPath,numProcs,thisVar)


fileName = 'output0.h5';
thisFile = [dataPath,fileName];
procID  = hdf5read(thisFile,'procID');
fileinfo = hdf5info(thisFile);
dsets=fileinfo.GroupHierarchy.Datasets;

nXg = 2;
try nXg = hdf5read(thisFile,'nXg');
end

%Name1 = fileinfo.GroupHierarchy.Datasets(1).Name

numVars = length(dsets);

Rank = -1;
for n=1:numVars
    thisVarName = dsets(n).Name(2:end);
    
    if(strcmp(thisVar,thisVarName))
        Rank = dsets(n).Rank;
      %  display(Rank);
        break;
    end
    
end
if(Rank==-1)
    disp([thisvar,' not found in output0.h5']);
end


%%%   get grid information from proc 0 output
%

if(Rank==0 || strcmp(thisVar,'Zcc') || strcmp(thisVar,'Zce') ...
           || strcmp(thisVar,'tout')) % scalar
    
   varData = hdf5read(thisFile,thisVar);
   
elseif(Rank==1) % grid vector
    
   varData = hdf5read(thisFile,thisVar);
   Xce = hdf5read(thisFile,'Xce');
   Xce2 = -100;
   try Xce2 = hdf5read(thisFile,'Xce2');
   end
   offset_ce = 0;
   if(length(varData)==length(Xce))
       offset_ce = 1;
   end
   if(length(varData)==length(Xce2))
       offset_ce = -1;
   end
   for j=1:numProcs-1
      fileName = ['output',num2str(j),'.h5'];
      thisFile = [dataPath,fileName];
      thisProcData = hdf5read(thisFile,thisVar);
    
      max0 = length(varData)-nXg+offset_ce;
      thisLength = length(thisProcData);
      varData(max0+1:max0+thisLength-nXg) = thisProcData(nXg+1:end);
      
      %max0 = length(varData)-2*nXg+offset_ce;
      %thisLength = length(thisProcData);
      %varData(max0+1:max0+thisLength) = thisProcData;
   end   
 
elseif(Rank==2) % vector Data growing in time
    
   varData = hdf5read(thisFile,thisVar);
   Xce = hdf5read(thisFile,'Xce');
   Xce2 = 100;
   try Xce2 = hdf5read(thisFile,'Xce2');
   end
   offset_ce = 0;
   if(length(varData(:,1))==length(Xce))
       offset_ce = 1;
   end
   if(length(varData(:,1))==length(Xce2))
       offset_ce = -1;
   end
   for j=1:numProcs-1
      fileName = ['output',num2str(j),'.h5'];
      thisFile = [dataPath,fileName];
      thisProcData = hdf5read(thisFile,thisVar);
    
      max0 = length(varData(:,1))-nXg+offset_ce;
      thisLength = length(thisProcData(:,1));
      %index0 = max0+1
      %index1 = max0-nXg+thisLength
      varData(max0+1:max0+thisLength-nXg,:) = thisProcData(nXg+1:end,:);
      
      %max0 = length(varData(:,1))-2*nXg+offset_ce;
      %thisLength = length(thisProcData(:,1));
      %varData(max0+1:max0+thisLength,:) = thisProcData;
   end
   
end
    
    
    
    
    
    



end