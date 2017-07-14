import CombineHarvester.CombineTools.ch as ch
import CombineHarvester.CombinePdfs.morphing as morphing
import ROOT
import glob,datetime,os,re

def harvest(subdir, mwrange, charge='both'):
    cmb = ch.CombineHarvester()
    
    # Read all the cards.
    # CH stores metadata about each object (Observation, Process, Systematic),
    # this is extracted from the card names with some regex
    for card in glob.glob(subdir+'/wenu_mass*.txt'): #.format(ch='' if charge == 'both' else 'plus*' if charge == 'plus' else 'minus*')):
        cmb.QuickParseDatacard(card, """wenu_mass(?<MASS>\d+)_$CHANNEL.card.txt""")
    
    # Need a unqiue bin name for each plus/minus,pt and eta combination
    # We extracted this part of the datacard name into the channel variable above,
    # so can just copy it and override the specific bin name that was in all the cards
    cmb.ForEachObj(lambda obj: obj.set_bin(obj.channel()))
 
    # We'll have three copies of the observation, one for each mass point.
    # Filter all but one copy.
    #cmb.FilterObs(lambda obj: obj.mass() != ('%d' % self.cmass))
    cmb.FilterObs(lambda obj: obj.mass() != '19')
    
    # Create workspace to hold the morphing pdfs and the mass
    w = ROOT.RooWorkspace('morph', 'morph')
    mass = w.factory('mw[{mwrange}]'.format(mwrange=mwrange))
    
    # BuildRooMorphing will dump a load of debug plots here
    debug = ROOT.TFile(subdir+'/debug.root', 'RECREATE')
    
    # Run for each bin,process combination (only for signal!)
    for b in cmb.bin_set():
        for p in cmb.cp().bin([b]).signals().process_set():
            morphing.BuildRooMorphing(w, cmb, b, p, mass, verbose=True, file=debug)
    
    # Just to be safe
    mass.setConstant(True)
    
    # Now the workspace is copied into the CH instance and the pdfs attached to the processes
    # (this relies on us knowing that BuildRooMorphing will name the pdfs in a particular way)
    cmb.AddWorkspace(w, True)
    cmb.cp().process(['W']).ExtractPdfs(cmb, 'morph', '$BIN_$PROCESS_morph', '')
    
    # Adjust the rateParams a bit - we currently have three for each bin (one for each mass),
    # but we only want one. Easiest to drop the existing ones completely and create new ones
    cmb.syst_type(['rateParam'], False)
    cmb.cp().process(['W']).AddSyst(cmb, 'norm_$BIN', 'rateParam', ch.SystMap()(1.00))
    
    # Have to set the range by hand
    for sys in cmb.cp().syst_type(['rateParam']).syst_name_set():
        cmb.GetParameter(sys).set_range(0.5, 1.5)
    
    # Print the contents of the model
    cmb.PrintAll()
    
    # Write out the cards, one per bin
    outdir=subdir+'/wenu_cards_morphed_{charge}'.format(charge=charge)
    writer = ch.CardWriter('$TAG/$BIN.txt', '$TAG/shapes.root')
    writer.SetVerbosity(1)
    writer.WriteCards(outdir, cmb)

date = datetime.date.today().isoformat()
date+='_charges'

card_dir = 'cards/cards_100717/'
subdirs = [x[0] for x in os.walk(card_dir)]

mwrange='0,38'
npoints = 39
central = 19

runHarvest = False
runBatch   = False
justHadd   = False
combineCards = False
runFit = True

input_dcs_alleta = ""
workspaces = []
for isub, subdir in enumerate(subdirs):
    if subdir == subdirs[0]: continue
    if 'wenu_cards_morphed' in subdir: continue
    name = subdir.split('/')[-1]
    if not 'eta_' in name: continue

    print '--------------------------------------------------------------------'
    print '- running for {mode} -----------------------------------------------'.format(mode=name)
    print '- in subdirectory {subdir} -----------------------------------------'.format(subdir=subdir)
    print '--------------------------------------------------------------------'
    #if name == 'full_3d': continue

    if runHarvest: 
        ## run the combine harvester which combines all the datacards etc.
        harvest(subdir,mwrange)

    target_dc = '{subdir}/wenu_cards_morphed_both/morphed_datacard_channel.txt'.format(subdir=subdir)
    target_ws = target_dc.replace('txt','root')
    workspaces.append(target_ws)

    if combineCards:
        ## running combineCards to make the combined plus+minus datacard
        if os.path.isfile(target_dc):
            print 'removing existing combined datacard first!'
            os.system('rm {dc}'.format(dc=target_dc) )
        dcs = os.listdir(subdir+"/wenu_cards_morphed_both/")
        input_dcs=" ".join(["%s=%s" % (os.path.splitext(dc)[0],subdir+"/wenu_cards_morphed_both/"+dc) for dc in dcs if "txt" in dc])
        input_dcs_alleta += " "+input_dcs

        print 'running combineCards.py'
        combineCardsCmd = 'combineCards.py {dcs} >& {target_dc}'.format(dcs=input_dcs, target_dc=target_dc)
        print combineCardsCmd
        ## run combineCards and make the workspace
        os.system(combineCardsCmd )
        print 'running text2workspace'
        t2wCmd = 'text2workspace.py {target_dc} '.format(subdir=subdir, target_dc=target_dc)
        print t2wCmd
        os.system(t2wCmd)

comb_dir = card_dir+'/comb'
if not os.path.exists(comb_dir): os.mkdir(comb_dir)
comb_dc = comb_dir+"/morphed_datacard_comb.txt"
comb_ws = comb_dc.replace('txt','root')
workspaces.append(comb_ws)
        
if combineCards:
    if os.path.isfile(comb_dc):
        print 'removing existing combined datacard first!'
        os.system('rm {dc}'.format(dc=comb_dc) )

    print 'running combineCards.py'
    combineCardsCmd = 'combineCards.py {dcs} >& {target_dc}'.format(dcs=input_dcs_alleta, target_dc=comb_dc)
    print combineCardsCmd
    ## run combineCards and make the workspace
    os.system(combineCardsCmd)
    print 'running text2workspace'
    os.system('text2workspace.py %s' % comb_dc)

if runFit:
    for m,ws in enumerate(workspaces):
        print "===> RUN FIT FOR WORKSPACE: ",ws
        name = re.search('\S+eta\_(\S+)\/wenu\S+',ws).group(1)
        if name==None: name="comb"
        ## constructing the command
        combine_base  = 'combine -t -1 -M MultiDimFit --setPhysicsModelParameters mw={central},r=1 --setPhysicsModelParameterRanges mw={mwrange} '.format(central=central,mwrange=mwrange)
        combine_base += ' --redefineSignalPOIs=mw --algo grid --points {npoints} {target_ws} '.format(npoints=npoints, target_ws=ws)
        
        saveNuisances = ''
        saveNuisances += ' --saveSpecifiedNuis {vs} '.format(vs=','.join('CMS_We_pdf'+str(i) for i in range(1,27)))
        
        run_combine_allUnc = combine_base + ' -n {date}_{name} {sn} '.format(date=date,name=name,sn=saveNuisances) 
        run_combine_noPdf  = combine_base + ' -n {date}_{name}_noPDFUncertainty --freezeNuisanceGroups pdfUncertainties '.format(date=date,name=name)
        run_combine_noPtW  = combine_base + ' -n {date}_{name}_noPTWUncertainty --freezeNuisances CMS_W_ptw '.format(date=date,name=name)
        run_combine_noEScale  = combine_base + ' -n {date}_{name}_noEScaleUncertainty --freezeNuisances CMS_We_elescale '.format(date=date,name=name)
        
        if runBatch:
            run_combine_allUnc += ' --job-mode lxbatch --split-points 10 --sub-opts="-q 8nh" --task-name {name}                  '.format(name=name)
            run_combine_noPdf  += ' --job-mode lxbatch --split-points 10 --sub-opts="-q 8nh" --task-name {name}_noPDFUncertainty '.format(name=name)
            run_combine_noPtW  += ' --job-mode lxbatch --split-points 10 --sub-opts="-q 8nh" --task-name {name}_noPtWUncertainty '.format(name=name)
            run_combine_noEScale  += ' --job-mode lxbatch --split-points 10 --sub-opts="-q 8nh" --task-name {name}_noEScaleUncertainty '.format(name=name)
            run_combine_allUnc  = 'combineTool.py ' + ' '.join(run_combine_allUnc.split()[1:])
            run_combine_noPdf   = 'combineTool.py ' + ' '.join(run_combine_noPdf .split()[1:])
            run_combine_noPtW   = 'combineTool.py ' + ' '.join(run_combine_noPtW .split()[1:])
            run_combine_noEScale   = 'combineTool.py ' + ' '.join(run_combine_noEScale .split()[1:])
        
        
        ## running combine once with the systematics and once without
        print '-- running combine command ------------------------------'
        print '---     with uncertainties: -----------------------------'
        print run_combine_allUnc        
        os.system(run_combine_allUnc)
        
        print '---     without PDF uncertainties: --------------------------'
        print run_combine_noPdf
        os.system(run_combine_noPdf )

        print '---     without PTW uncertainties: --------------------------'
        print run_combine_noPtW
        os.system(run_combine_noPtW )

        print '---     without electron energy scale uncertainties: --------------------------'
        print run_combine_noEScale
        os.system(run_combine_noEScale )
        
        impactBase = 'combineTool.py -M Impacts -n {date}_eta_{name} -d {target_ws} -m {mass}'.format(mass=m,date=date,name=name, target_ws=ws)
        impactBase += ' --setPhysicsModelParameters mw={central},r=1  --redefineSignalPOIs=mw --setPhysicsModelParameterRanges mw={mwrange} -t -1 '.format(central=central,mwrange=mwrange)
        impactInitial = impactBase+'  --robustFit 1 --doInitialFit '
        impactFits    = impactBase+'  --robustFit 1 --doFits '
        impactJSON    = impactBase+'  -o impacts_eta_{name}.json '.format(name=name)
        impactPlot    = 'plotImpacts.py -i impacts_eta_{name}.json -o impacts_eta_{name} --transparent'.format(name=name)
        
        # os.system(impactInitial)
        # os.system(impactFits   )
        # os.system(impactJSON   )
        # os.system(impactPlot   )



