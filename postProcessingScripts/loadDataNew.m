function[varData]=loadDataNew(dataPath,thisVar)


fileName = 'output0.h5';
thisFile = [dataPath,fileName];
procID  = hdf5read(thisFile,'procID');
fileinfo = hdf5info(thisFile);
dsets=fileinfo.GroupHierarchy.Datasets;


%%%   get the number of processors
%
fileNames = dir([dataPath,'output*.h5']);
numProcs = length(fileNames);


%%%   get the data size
%
nXg = 2;
try nXg = hdf5read(thisFile,'nXg');
end
nZg = 2;
try nZg = hdf5read(thisFile,'nZg');
end   


%%%   get X sizes per proc and total
%
Xcc = hdf5read(thisFile,'Xcc');
Nxsub_cc = length(Xcc);
Nxsub_ce = Nxsub_cc-1;
Nxsub_ce2 = Nxsub_cc+1;
%
Nx_cc = numProcs*(Nxsub_cc-2*nXg) + 2*nXg;
Nx_ce = Nx_cc-1;
Nx_ce2 = Nx_cc+1;

varData0 = hdf5read(thisFile,thisVar);

%Name1 = fileinfo.GroupHierarchy.Datasets(1).Name

numVars = length(dsets);

Rank = -1;
for n=1:numVars
    thisVarName = dsets(n).Name(2:end);
    
    if(strcmp(thisVar,thisVarName))
        Rank = dsets(n).Rank;
       % display(Rank);
        break;
    end
    
end
if(Rank==-1)
    disp([thisvar,' not found in output0.h5']);
end


%%%   get grid information from proc 0 output
%

if(Rank==0 || strcmp(thisVar,'Zcc') || strcmp(thisVar,'Zce') ...
           || strcmp(thisVar,'Zce2') || strcmp(thisVar,'tout')) % scalar
    
   varData = varData0; %hdf5read(thisFile,thisVar);
   
elseif(Rank==1) % X-grid vector
    
   varData = hdf5read(thisFile,thisVar);

   offset_ce = 0;
   Nxvar = Nx_cc;
   if(length(varData)==Nxsub_ce )
       offset_ce = 1;
       Nxvar = Nx_ce;
       Nxsub = Nxsub_ce;
   end
   if(length(varData)==Nxsub_ce2 )
       offset_ce = -1;
       Nxvar = Nx_ce2;
       Nxsub = Nxsub_ce2;
   end
   for j=1:numProcs-1
      fileName = ['output',num2str(j),'.h5'];
      thisFile = [dataPath,fileName];
      thisProcData = hdf5read(thisFile,thisVar);
      
      max0 = length(varData)-nXg+offset_ce;
      thisLength = length(thisProcData);
      varData(max0+1:max0+thisLength-nXg) = thisProcData(nXg+1:end);
   end
   
elseif(Rank==2) % matrix Data not growing in time
    
   offset_ce = 0;
   Nxvar = Nx_cc;
   Nxsub = Nxsub_cc;
   if(length(varData0(1,:))==Nxsub_ce )
       offset_ce = 1;
       Nxvar = Nx_ce;
       Nxsub = Nxsub_ce;
   end
   if(length(varData0(1,:))==Nxsub_ce2 )
       offset_ce = -1;
       Nxvar = Nx_ce2;
       Nxsub = Nxsub_ce2;
   end
   
   varData = zeros(length(varData0(:,1)),Nxvar);
   varData(:,1:length(varData0(1,:))) = varData0;
   
   i0 = length(varData0(1,:))-nXg+1+offset_ce;
   for j=1:numProcs-1
      fileName = ['output',num2str(j),'.h5'];
      thisFile = [dataPath,fileName];
      thisProcData = hdf5read(thisFile,thisVar);
      
      i1 = i0 - 1 + Nxsub - nXg;
      varData(:,i0:i1) = thisProcData(:,nXg+1:end);
      i0 = i1 - nXg + 1 + offset_ce;
   end
 
elseif(Rank==3) % matrix Data growing in time
    
   offset_ce = 0;
   Nxvar = Nx_cc;
   Nxsub = Nxsub_cc;
   if(length(varData0(1,:,1))==Nxsub_ce )
       offset_ce = 1;
       Nxvar = Nx_ce;
       Nxsub = Nxsub_ce;
   end
   if(length(varData0(1,:,1))==Nxsub_ce2 )
       offset_ce = -1;
       Nxvar = Nx_ce2;
       Nxsub = Nxsub_ce2;
   end
   
   varData = zeros(length(varData0(:,1,1)),Nxvar,length(varData0(1,1,:)));
   varData(:,1:length(varData0(1,:,1)),:) = varData0;
   
   i0 = length(varData0(1,:,1))-nXg+1+offset_ce;
   for j=1:numProcs-1
      fileName = ['output',num2str(j),'.h5'];
      thisFile = [dataPath,fileName];
      thisProcData = hdf5read(thisFile,thisVar);
      
      i1 = i0 - 1 + Nxsub - nXg;
      varData(:,i0:i1,:) = thisProcData(:,nXg+1:end,:);
      i0 = i1 - nXg + 1 + offset_ce;
   end
   
end
    
    
    

end
