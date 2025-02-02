clc;close all; clear all
%Use these Keywords to find related article to your topics (This is case
%sensetive! so make sure you 
KeyWords=["Retina";"retina";];
KeyWords=[KeyWords;upper(KeyWords)];
clc;ArticleList=strings(500,1);counter=1;


%% Nature neuroscience
clc;
url = "https://www.nature.com/search?q=neuroscience&journal=ncomms&order=date_desc";
code = webread(url);


for p = 1 : length(KeyWords)
    %k will hold indexes of key words for the serch:
    if contains(code,KeyWords(p))
        k = strfind(code,(KeyWords(p)));
        
        for kk=1:length(k)
            flag= false;
            %Finiding out the title of the article (By general - site specific cheractaristics:
            
            ind1=strfind(code,'<a href=');
            temp=ind1-k(kk);
            temp=find(temp<0);
            ind1=ind1(temp(end));

            ind2=strfind(code,'class="');
            temp=ind2-k(kk);
            temp=find(temp<0);
            ind2=ind2(temp(end-1));
            
            ArticleList(counter)=['https://www.nature.com',(code(ind1+9:ind2-1))];
            counter=counter+1;
        end
    end
end

%% Cell press
url = "https://www.cell.com/neuron/newarticles";
code = webread(url);


for p = 1 : length(KeyWords)
    %k will hold indexes of key words for the serch:
    if contains(code,KeyWords(p))
        k = strfind(code,(KeyWords(p)));
        
        for kk=1:length(k)
            flag= false;
            %Finiding out the title of the article (By general - site specific cheractaristics:
            
            ind1=strfind(code,'<a href="');
            temp=ind1-k(kk);
            temp=find(temp<0);
            ind1=ind1(temp(end));

            ind2=strfind(code,'" title=');
            temp=ind2-k(kk);
            temp=find(temp<0);
            ind2=ind2(temp(end));
            
            ArticleList(counter)=['https://www.cell.com/',(code(ind1+9:ind2-1))];
            counter=counter+1;
        end
    end
end

%% PNAS
url = "https://www.pnas.org/topic/neuro";
code = webread(url);


for p = 1 : length(KeyWords)
    %k will hold indexes of key words for the serch:
    if contains(code,KeyWords(p))
        k = strfind(code,(KeyWords(p)));
        
        for kk=1:length(k)
            flag= false;
            %Finiding out the title of the article (By general - site specific cheractaristics:
            
            ind1=strfind(code,'<a href="');
            temp=ind1-k(kk);
            temp=find(temp<0);
            ind1=ind1(temp(end));
      
            ind2=strfind(code,'" title=');
            temp=ind2-k(kk);
            temp=find(temp<0);
            ind2=ind2(temp(end));
            
            ArticleList(counter)=['https://www.pnas.org',(code(ind1+9:ind2-1))];
            counter=counter+1;
        end
    end
end

%% eLife
url = "https://elifesciences.org/subjects/neuroscience";
code = webread(url);


for p = 1 : length(KeyWords)
    %k will hold indexes of key words for the serch:
    if contains(code,KeyWords(p))
        k = strfind(code,(KeyWords(p)));
        
        for kk=1:length(k)
            flag= false;
            %Finiding out the title of the article (By general - site specific cheractaristics:
            
            ind1=strfind(code,'<a href="');
            temp=ind1-k(kk);
            temp=find(temp<0);
            ind1=ind1(temp(end));

            ind2=strfind(code,'"   class="');
            temp=ind2-k(kk);
            temp=find(temp<0);
            ind2=ind2(temp(end));
            
            ArticleList(counter)=['https://elifesciences.org',(code(ind1+9:ind2-1))];
            counter=counter+1;
        end
    end
end

%% jneurosci
url = "https://www.jneurosci.org/content/early/recent";
code = webread(url);


for p = 1 : length(KeyWords)
    %k will hold indexes of key words for the serch:
    if contains(code,KeyWords(p))
        k = strfind(code,(KeyWords(p)));
        
        for kk=1:length(k)
            flag= false;
            %Finiding out the title of the article (By general - site specific cheractaristics:
            
            ind1=strfind(code,'<a href="');
            temp=ind1-k(kk);
            temp=find(temp>1);
            ind1=ind1(temp(1));
            
            ind2=strfind(code,'" class=');
            temp=ind2-k(kk);
            temp=find(temp>1);
            ind2=ind2(temp(1));
            
            ArticleList(counter)=['https://www.jneurosci.org',(code(ind1+9:ind2-1))];
            counter=counter+1;
        end
    end
end

%% BioArchive
HighURL=["https://www.biorxiv.org/collection/neuroscience";...
    "https://www.biorxiv.org/collection/neuroscience?page=1";...
    "https://www.biorxiv.org/collection/neuroscience?page=3";...
    "https://www.biorxiv.org/collection/neuroscience?page=4"];

for q=1:length(HighURL)
    url = HighURL(q);
    code = webread(url);
    
    
    for p = 1 : length(KeyWords)
        %k will hold indexes of key words for the serch:
        if contains(code,KeyWords(p))
            k = strfind(code,(KeyWords(p)));
            
            for kk=1:length(k)
                flag= false;
                %Finiding out the title of the article (By general - site specific cheractaristics:
                
                ind1=strfind(code,'<a href="');
                temp=ind1-k(kk);
                temp=find(temp<0);
                ind1=ind1(temp(end));
                
                ind2=strfind(code,'" class="');
                temp=ind2-k(kk);
                temp=find(temp<0);
                ind2=ind2(temp(end));
                
                ArticleList(counter)=['https://www.biorxiv.org',(code(ind1+9:ind2-1))];
                counter=counter+1;
            end
        end
    end
end
%%
ArticleList=unique(ArticleList);
 filename = 'New Article!!.xlsx';
 writematrix(ArticleList,filename)
 
%  %%
%   
%  clc;close all; clear all
% KeyWords=["dopamine";"striatum";"striatal";"da ";"SNc ";"dRN ";"5-HT";...
%     "basal ganglia";"pacap";"synaptic";"muscarin"];
% KeyWords=[KeyWords;upper(KeyWords)];
% clc;ArticleList=strings(500,1);counter=1;
% 
% 
% %% Science
% url = "https://www.science.org/action/doSearch?AllField=neuroscience&SeriesKey=science&startPage=0&sortBy=Earliest";
% code = webread(url);
% 
% 
% for p = 1 : length(KeyWords)
%     %k will hold indexes of key words for the serch:
%     if contains(code,KeyWords(p))
%         k = strfind(code,(KeyWords(p)));
%         
%         for kk=1:length(k)
%             flag= false;
%             %Finiding out the title of the article (By general - site specific cheractaristics:
%             
%             ind1=strfind(code,'<a href="');
%             temp=ind1-k(kk);
%             temp=find(temp<0);
%             ind1=ind1(temp(end));
%             
%             %Original
% %             ind1=strfind(code,'<a href=');
% %             temp=ind1-k(kk);
% %             ind2=find(temp<0);
% %             ind3=ind1(ind2(end));
%             
%             ind2=strfind(code,'" class="');
%             temp=ind2-k(kk);
%             temp=find(temp<0);
%             ind2=ind2(temp(end));
%             
%             ArticleList(counter)=['https://www.science.org/',(code(ind1+9:ind2-1))];
%             counter=counter+1;
%         end
%     end
% end
% 
% ArticleList=unique(ArticleList);
% %%
%  filename = 'testdata.xlsx';
%  writematrix(ArticleList,filename)
 
 