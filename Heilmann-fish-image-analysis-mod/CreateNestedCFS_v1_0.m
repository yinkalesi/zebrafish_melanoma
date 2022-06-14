%% CreateNestedCFS_v1_0
%  Version 1.0
%  Authors: Adeyinka Lesi
%  Date: 2/13/20
%% Version History
%  1.0: to resolve memory handling issues, will create a nested version of
%  CreateFishStructFunction


template = 'CreateFishStructFunction_v4_1';
output = 'CreateFishStructNested_v4_1';

% 1) find functions required
% search file for functions using BYTIME
search_text = '([%]?)(.*) [=]? (\w+)\((.*)BYTIME[,]?([^\(\)]*)[\)\.]';

fid = fopen([template '.m']);
C = textscan(fid,'%s','Delimiter','\n');
fclose(fid);

% extract file names and locations of those functions
[~,~,~,~,tokens]=regexp(C{1}(:),search_text,'once');
match_line = find(~cellfun('isempty',tokens));

% get file_names of functions from tokens var
function_files = cell(1,length(match_line));
for f = 1:length(match_line)
    function_files{f} = tokens{match_line(f)}{3};
end

% 2) find insertion point in template
% find end statement if present (take care due to loops and if statements!)
last_text_line = find(~cellfun('isempty',regexp(C{1}(:),'.+','once')),1,'last');
if(strcmp(C{1}(last_text_line),'end'))
    end_line = last_text_line;
else
    end_line = length(C{1})+1;
end

% 3) construct file
% open new file
fid = fopen([output '.m'],'w');

% put code in, make sure to rewrite the function calls so BYTIME is not
% passed
func_line = find(~cellfun('isempty',regexp(C{1}(:),'^\W*function','once')),1,'first');
for line = 1:end_line-1
    if line==func_line
        new_func_line = regexprep(C{1}{func_line},'= \w+\(',['= ' output '(']);
        fprintf(fid,'%% ***MODIFIED using %s***\n',mfilename);
        fprintf(fid,'%s\n',new_func_line);
    elseif ~any(line==match_line)
        % just copy and past
        fprintf(fid,'%s\n',C{1}{line});
    else
        fprintf(fid,'%% ***MODIFIED using %s***\n',mfilename);
        [~,~,~,~,~,~,split] = regexp(tokens{line}{2},'BYTIME[,]?','once');
        output_text = join(split);
        if(regexp(output_text,'\w+','once'))
            fprintf(fid,'%s %s = %s(%s%s);\n',tokens{line}{1},output_text,tokens{line}{3},tokens{line}{4},tokens{line}{5});
        else
            fprintf(fid,'%s %s(%s%s);\n',tokens{line}{1},tokens{line}{3},tokens{line}{4},tokens{line}{5});
        end
    end
end

% write in nested functions, making sure BYTIME is never an input or output
for f = 1:length(function_files)
    fid2 = fopen([function_files{f} '.m']);
    C2 = textscan(fid2,'%s','Delimiter','\n');
    fclose(fid2);
    
    % find the function definition
    func_line = find(~cellfun('isempty',regexp(C2{1}(:),'^\W*function','once')),1,'first');
    if(isempty(func_line))
        error('Expected %s to contain a function\n',function_files{f});
    end
    % reconstruct without BYTIME if necessary
    [~,~,~,~,toks2] = regexp(C2{1}{func_line},['\W*function (.*) [=]? ' function_files{f} '\((.*)BYTIME[,]?(.*)[\)\.]'],'once');
    
    if(isempty(toks2))
        new_func_text = C2{1}{func_line};
    else
        [~,~,~,~,~,~,split] = regexp(toks2{1},'BYTIME[,]?','once');
        new_func_text = sprintf('function %s = %s(%s%s)',join(split),function_files{f},toks2{2},toks2{3});
    end
    
    fprintf(fid,'\n\n%% ***INSERTED using %s***\n',mfilename);
    for f2 = 1:func_line-1
        fprintf(fid,'%s\n',C2{1}{f2});
    end
    fprintf(fid,'%s\n',new_func_text);
    for f2 = func_line+1:length(C2{1})
        fprintf(fid,'%s\n',C2{1}{f2});
    end
    % make sure last text is an end statement
    last_text2 = find(~cellfun('isempty',regexp(C2{1}(:),'.+','once')),1,'last');
    if(~strcmp(C2{1}(last_text2),'end'))
        fprintf(fid,'\nend\n');
    end
end
    

% finish up
fprintf(fid,'\n\nend');
fclose(fid);