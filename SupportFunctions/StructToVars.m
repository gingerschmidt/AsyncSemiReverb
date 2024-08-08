function StructToVars(varsStruct)
    
  % Unpacks a struct containing variables. 
  % 
  % Author: Néstor Uribe-Patarroyo (1) 
  % 1. Wellman Center for Photomedicine, Harvard Medical School, Massachusetts
  % General Hospital, 40 Blossom Street, Boston, MA, USA
  % <uribepatarroyo.nestor@mgh.harvard.edu>
  
  varNames = fieldnames(varsStruct);    
  for k = 1:numel(varNames)
    assignin('caller', varNames{k}, varsStruct.(varNames{k}));
  end
  
end

