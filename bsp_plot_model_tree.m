function varargout=bsp_plot_model_tree(what,varargin)

% usage: bsp_plot_model_tree('plot_GLM_Physio_full_model')

subj_name = {'S98','S97','S96','S95','S01','S03','S04','S07'};

switch(what)
    
case 'plot_GLM_Physio_full_model'
    roi = {'simulate'};
    what = 'R'; % what to plot - here correlation on
    sn = [1:8];
    vararginoptions(varargin,{'roi','what','sn'});

    D=load('test_GLM_physio_simulate_tikhonov.mat');

    num_subj = length(sn);
    for s=1:num_subj
        figure;
        for r=1:length(roi)
            subplot(1,2,r);
            
            tmp = D.(what)(find(D.sn==sn(s)));
            indx = find(contains(D.roi(find(D.sn==sn(s))),roi{r}));
            avg = splitapply(@mean,tmp(indx),D.model(1:numel(indx)));
            
            start = [1 1 1 1 1 2 2 2 2 3 3 3 3 4 4 4 4 5 5 5 5 6 6 6 6 7 7 7 8 8 8 9 9 9 10 10 10 11 11 11 12 12 12 13 13 13 14 14 14 15 15 15 16 16 16 ...
                17 17 18 18 19 19 20 20 21 21 22 22 23 23 24 24 25 25 26 26 27 28 29 30 31];
            terminate = [2 3 4 5 6 7 8 9 10 7 11 12 13 8 11 14 15 9 12 15 16 10 13 15 16 17 18 19 17 20 21 18 20 22 19 21 22 17 23 24 18 23 25 19 24 25 18 23 26 21 24 26 22 25 26 ...
                27 28 27 29 28 29 27 30 28 30 29 30 27 31 28 31 29 31 30 31 32 32 32 32 32];
            
            avgTemp = [0; avg];
            
            for i=1:length(start)
                weights(i) = avgTemp(terminate(i))-avgTemp(start(i));
            end
            code = round(weights,3);
            
            EdgeTable = table([start' terminate'],weights', code','VariableNames',{'EndNodes','Weight','Code'});
            
            names = {'Null','Task','RetHR','RetRESP','HR','RV','Task+RetHR','Task+RetRESP','Task+HR','Task+RV','RetHR+RetRESP','RetHR+HR','RetHR+RV','RetRESP+HR','RetRESP+RV','HR+RV', ...
                    'Task+RetHR+RetRESP','Task+RetHR+HR','Task+RetHR+RV','Task+RetRESP+HR','Task+RetRESP+RV','Task+HR+RV','RetHR+RetRESP+HR','RetHR+RetRESP+RV', ...
                    'RetHR+HR+RV','RetRESP+HR+RV','Task+RetHR+RetRESP+HR','Task+RetHR+RetRESP+RV','Task+RetHR+HR+RV','Task+RetRESP+HR+RV','RetHR+RetRESP+HR+RV','All'}';
            namesCode = [1:1:32]';
            nodeWeights = [0; avg];
            NodeTable = table(names,namesCode,nodeWeights,'VariableNames',{'Name','NameCode','NodeWeight'});
            
            G = digraph(EdgeTable,NodeTable);
            
            G.Nodes.NodeWeight(isnan(G.Nodes.NodeWeight)) = 0;
            G.Edges.Weight(isnan(G.Edges.Weight)) = 1e-6;
            
            G.Nodes.NodeColors = G.Nodes.NodeWeight;
            
            xData = [7,4,5.5,7,8.5,10,0.25,1.75,3.25,4.75,6.25,7.75,9.25,10.75,12.25,13.75,0.25,1.75,3.25,4.75,6.25,7.75,9.25,10.75,12.25,13.75,4,5.5,7,8.5,10,7];
            yData = [6,5,5,5,5,5,4,4,4,4,4,4,4,4,4,4,3,3,3,3,3,3,3,3,3,3,2,2,2,2,2,1];
            %figure; g = plot(G,'NodeLabel',G.Nodes.NameCode,'EdgeLabel',G.Edges.Code,'Layout','layered');
            g = plot(G,'NodeLabel',G.Nodes.NameCode,'XData',xData','YData',4*yData);
            title(sprintf('%s %s %s',subj_name{sn(s)},roi{r},D.method{1}));
            g.NodeCData = G.Nodes.NodeColors;
            g.MarkerSize = 10;
            g.LineWidth = abs(100*G.Edges.Weight);
            g.EdgeColor = [0.5 0.5 0.5];
            colorbar
        end
    end
end

