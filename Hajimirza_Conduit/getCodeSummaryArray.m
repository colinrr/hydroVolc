function plotCodes = getCodeSummaryArray(outcomeCodes)
%  plotCodes = getCodeSummaryArray(outcomeCodes)
% temporary function to retrieve most important code results
%
%   Problem: outcomeCodes.maxR is spitting out zero codes instead of
%   negatives when runs fail, but the correct codes are recorded in
%   outcomeCodes.all. So this function searches out such cases and
%   populates the correct codes into an array as a temporary fix.

    knownFailsToMap = {'Check Pg, T'}; % update code assignment below if needed

    maxCode = cellfun(@max,outcomeCodes.all);
    % Account for new to-be-mapped errors
    unmappedFails = (maxCode == -22);
    umi = find(unmappedFails);
    for um = 1:sum(unmappedFails(:))
        if ismember(outcomeCodes.failMsg(umi(um)).ME.message,knownFailsToMap)
            maxCode(umi(um)) = -21; % update as needed if knownFailsToMap is updated
        end
    end
    % Update
    codeZero = or(outcomeCodes.maxR < 1,isnan(outcomeCodes.maxR));
    plotCodes = outcomeCodes.maxR;
    plotCodes(codeZero) = maxCode(codeZero);
    
end