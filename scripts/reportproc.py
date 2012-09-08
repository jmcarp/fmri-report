
# Import Python libraries
import re
import warnings
import operator
import collections
from pylab import *
import scipy.stats

# Set up rpy2
import rpy2.robjects as ro
from rpy2.robjects.packages import importr
gr = importr('grDevices')

# Paths
homedir = '/Users/jmcarp/Dropbox/projects/fmri-report'
datadir = '%s/data' % (homedir)
scriptdir = '%s/scripts' % (homedir)
figsdir = '%s/figs' % (homedir)

# Import R libraries
ro.r('library(RColorBrewer)')
ro.r('library(plotrix)')
ro.r('library(pwr)')
ro.r('source("%s/modbar.R")' % (scriptdir))

# 
esizes = [0.8, 1.2, 1.6]
elabels = ['expression(paste(italic("d"), "=%0.1f", sep=""))' \
  % esize for esize in esizes]
estrs = [str(esize).replace('.', '') for esize in esizes]

# Plotting defaults
cex_lab = 1.75
cex_axis = 1.5
cex_leg = 1.75
cex_main = 2.25
mar = [5.6, 5.1, 4.1, 2.1]
mgpx = [4.25, 1, 0]
mgpy = [3, 1, 0]

tsvfile = '%s/fmri-report.tsv' % (datadir)

mandsteps = [ 
  'edes', 'subjs', 'scan', 'func', 
  'smod', 'gmod', 'misc'
]
boolsteps = [ 
  'eventopt', 'power', 'anat', 'coplanar', 
  'slicetime', 'mc', 'coreg', 'skullstrip', 
  'segment', 'norm', 'smooth', 'detrend', 
  'unwarp', 'scale', 'filter', 'motreg', 
  'acf', 'mccorrect', 'roi', 'fig', 'tab'
]
allsteps = mandsteps + boolsteps

legpos = {
  'mod' : 'topleft',
  'proc' : 'topleft',
  'displ' : 'topleft'
}

ymax = {
  'misc-softcom' : 80
}

merge = {
  'proc-norm-tgtvol' : {
    'idealab t1 template' : 'custom template',
    'average t2' : 'custom template',
    'custom atlas' : 'custom template',
    'custom epi atlas' : 'custom template',
    'mni305 atlas' : 'mni305 template',
    'pals atlas' : 'pals template',
    'mni152 atlas' : 'mni152 template',
    'mni251 t1 template' : 'mni152 template'
  },
  'misc-anatloc' : {
    'talairach daemon' : 'talairach atlas',
    'anatomy' : 'anatomical landmarks',
    'tracing' : 'anatomical landmarks',
    'sulcal landmarks' : 'anatomical landmarks'
  },
  'mod-gmod-mccorrect-method' : {
    'alphasim' : 'monte carlo'
  },
  'mod-smod-filter-filttype' : {
    'bandpass' : 'band-pass',
    'high-passs' : 'high-pass'
  },
  'misc-jtitle' : {
    'proceedings of the national academy of sciences of the united states of america' : 'pnas',
  }
}

stepgrp = {
  'des' : [ 'edes', 'subjs' ],
  'acq' : [ 'func', 'anat', 'coplanar' ],
  'proc' : [ 'slicetime', 'mc', 'coreg', 'norm', 'smooth' ],
  'mod' : [ 'smod', 'gmod', 'roi' ],
  'displ' : [ 'fig', 'tab' ]
}

stepname = {
  'edes' : 'Task design',
  'subjs' : 'Human subjects',
  'func' : 'Functional scans',
  'anat' : 'Anatomical scans',
  'coplanar' : 'Auxiliary scans',
  'slicetime' : 'Slice-timing correction',
  'mc' : 'Motion correction',
  'coreg' : 'Coregistration',
  'norm' : 'Spatial normalization',
  'smooth' : 'Spatial smoothing',
  'smod' : 'Within-subjects modeling',
  'gmod' : 'Between-subjects modeling',
  'roi' : 'Region of interest analysis',
  'fig' : 'Figures',
  'tab' : 'Tables'
}

grpname = {
  'des' : 'Design',
  'acq' : 'Acquisition',
  'proc' : 'Pre-processing',
  'mod' : 'Modeling',
  'displ' : 'Display'
}

def ezread():
  
  rep = readrep(tsvfile)
  rep = augrep(rep)
  rep = mergerep(rep)

  return rep

def flexplots():
  
  rep = readrep(tsvfile)
  rep = augrep(rep)
  rep = mergerep(rep)

  plotsmooth(rep, '%s/flex-smooth.pdf' % (figsdir))
  plothpf(rep, cutoff=500, binsize=20, outname='%s/flex-hpf.pdf' % (figsdir))
  plotext(rep, exttype='vol', outname='%s/flex-mccorrect-ext.pdf' % (figsdir))

  plotopts(rep, 'norm', 'proc-norm-tgtvol', minct=2, 
    outname='%s/flex-norm-tgtvol.pdf' % (figsdir))
  plotopts(rep, 'gmod', 'mod-gmod-mccorrect-height', minct=2,
    main='expression(bold("Inference: Height"))',
    outname='%s/flex-mccorrect-height.pdf' % (figsdir))
  plotopts(rep, 'misc', 'misc-softpck', minct=2, 
    main='expression(bold("Software packages"))',
    outname='%s/flex-misc-softpck.pdf' % (figsdir))
  plotopts(rep, 'misc', 'misc-softcom', minct=2, 
    main='expression(bold("Software versions"))',
    outname='%s/flex-misc-softcom.pdf' % (figsdir))

def makeplots():
  
  rep = readrep(tsvfile)
  rep = augrep(rep)

  boolinfo = procbool(rep)
  stepinfo = procsteps(rep)

  plotbool(boolinfo, rep, 'prop', 'TRUE', horiz=False,
    outname='%s/boolfig.pdf' % (figsdir)
  )

  for grp in stepgrp:
    plotmissing(stepinfo, rep, True, 'prop', 
      stepgrp[grp], main='expression(bold("%s"))' % grpname[grp],
      outname='%s/miss_%s.pdf' % (figsdir, grp))

def sampleplots():
  
  rep = readrep(tsvfile)
  rep = augrep(rep)
  subjinfo(rep)

  plotsubjs(rep, 100, binsize=4, 
    outname='%s/nsubj.pdf' % (figsdir))
  
  cumpwr1 = []
  cumpwr2 = []

  for eidx in range(len(esizes)):
    bins, cp1, cp2 = getcumpwr(rep, esizes[eidx], 0.001)
    cumpwr1.append(cp1)
    cumpwr2.append(cp2)
    #plotpwr(rep, esize, 0.001, binsize=0.05, 
    #  outname='%s/pwr-%s-001.pdf' % (figsdir, estrs[eidx]))

  plotcumpwr(bins, cumpwr1, 
    main='expression(bold("Power: One-group studies"))',
    outname='%s/pwr1.pdf' % (figsdir))
  plotcumpwr(bins, cumpwr2, 
    main='expression(bold("Power: Two-group studies"))',
    outname='%s/pwr2.pdf' % (figsdir))
  
def plotcumpwr(pwrbins, cumpwr, main=None, outname=None):
  
  # Open PDF
  if outname:
    gr.pdf(outname)
  
  # Set margins
  ro.r.par(mar=ro.FloatVector(mar))

  plotopts = {
    'cex.axis' : cex_axis,
    'cex.lab' : cex_lab,
  }

  ro.r.plot(1, 
    xlim=ro.r.c(0, 100), 
    ylim=ro.r.c(0, 1), 
    type='n',
    xlab='Power (%)',
    ylab='Proportion',
    **plotopts
  )
  
  cols = ['red', 'green', 'blue']

  for eidx in range(len(esizes)):
    ro.r.lines(
      ro.FloatVector(pwrbins),
      ro.FloatVector(cumpwr[eidx]),
      type='s',
      lwd=3,
      col=cols[eidx],
    )
  ro.r.abline(v=80, lty=2, lwd=3)
  
  if main:
    ro.r.title(main=ro.r(main), **{'cex.main' : cex_main})

  ro.r.legend('topright', ro.r('c(%s)' % ', '.join(elabels)),
    fill=ro.StrVector(cols), cex=cex_leg, inset=0.025, bg='white')

  # Close PDF
  if outname:
    gr.dev_off()


def pwr80(esize0, esizeinc, sig):
  
  rep = readrep(tsvfile)
  rep = augrep(rep)
  subjinfo(rep)
  
  p8res = {}
  
  pwr = 0
  esize = esize0
  while pwr < 0.8:
    esize += esizeinc
    getpwr(rep, esize, sig)
    pwrtmp = [r['subjinfo']['pwr'] for r in rep
      if 'subjinfo' in r and r['subjinfo']['gt1'] 
      and r['subjinfo']['ngrp'] == 1]
    pwr = float(len([p for p in pwrtmp if p > 0.8])) / len(pwrtmp)
  
  p8res['g1'] = esize

  pwr = 0
  esize = esize0
  while pwr < 0.8:
    esize += esizeinc
    getpwr(rep, esize, sig)
    pwrtmp = [r['subjinfo']['pwr'] for r in rep
      if 'subjinfo' in r and r['subjinfo']['gt1'] 
      and r['subjinfo']['ngrp'] == 2]
    pwr = float(len([p for p in pwrtmp if p > 0.8])) / len(pwrtmp)

  p8res['g2'] = esize

  return p8res

divcols = [ 'BrBG', 'PuOr', 'PiYG' ]

grpcol = collections.OrderedDict([
  ('des', 'Blues'), 
  ('acq', 'Oranges'),
  ('proc', 'Purples'),
  ('mod', 'Greens'),
  ('displ', 'PiYG')
])

def getspecmid(specname):
  if specname in divcols:
    cols = ro.r('brewer.pal(10, "%s")' % (specname))
  else:
    cols = ro.r('brewer.pal(5, "%s")' % (specname))
  return cols[2]

grppal = collections.OrderedDict([
  (key, getspecmid(grpcol[key])) for key in grpcol
])

def plotopts(rep, step, fld, usestrict=True, minct=5, rules=[], main=None, outname=None, taboutname=None):
  
  increp = [r for r in rep 
    if False not in [rule(r) for rule in rules]
  ]
  stepinfo = procsteps(increp)

  if 'strict' in stepinfo[step][fld]:
    if usestrict:
      rule = 'strict'
    else:
      rule = 'lenient'
  else:
    rule = 'mand'

  fldinfo = stepinfo[step][fld][rule]
  incopts = [opt for opt in fldinfo['count']
    if opt not in ['missing', 'n/a']]
  fldopts = [opt for opt in incopts if
    fldinfo['count'][opt] >= minct]
  othopts = [opt for opt in incopts if
    fldinfo['count'][opt] < minct]
  
  grays = ro.r('gray.colors(2)')
  cols = [grays[1]] * len(fldopts)

  fldvals = [fldinfo['count'][opt] for opt in fldopts]
  if len(othopts) > 0:
    othct = sum([fldinfo['count'][opt] for opt in othopts])
    fldopts.append('other')
    fldvals.append(othct)
    cols.append(grays[0])
  
  baropts = {
    'cex.axis' : cex_axis,
    'cex.lab' : cex_lab,
    'ylab' : 'Count'
  }
  ro.r.barplot(ro.IntVector(fldvals), col=ro.StrVector(cols),
    **baropts)
  
  ticks = ro.r.axTicks(2)
  ymax = max(fldvals)
  if ticks[-1] < ymax:
    ymax = ticks[-1] + (ticks[-1] - ticks[-2])
    ylim = 'c(%f, %f)' % (ymax * -0.01, ymax)
  else:
    ylim = 'NULL'
  
  # Open PDF
  if outname:
    gr.pdf(outname)
  
  # Set margins
  ro.r.par(mar=ro.FloatVector(mar))
  
  baropts = {
    'cex.axis' : cex_axis,
    'cex.lab' : cex_lab,
    'ylab' : 'Count',
    'ylim' : ro.r(ylim),
  }
  mid = ro.r.barplot(ro.IntVector(fldvals), col=ro.StrVector(cols),
    **baropts)
  
  for optidx in range(len(fldopts)):
    src = re.search('(\w+) ([\<\>]) (.*?$)', fldopts[optidx])
    if src:
      fldopts[optidx] = re.sub('(\w+) ([\<\>]) (.*?$)', 'paste(italic("\\1"), "\\2", "\\3")', fldopts[optidx])
    else:
      fldopts[optidx] = '"%s"' % fldopts[optidx]
  a2z = [chr(x) for x in range(ord('a'), ord('z') + 1)]
  labstr = 'c(' + \
    ', '.join(['expression(phantom("%s")*%s)' % (a2z, opt) for opt in fldopts]) + \
    ')'
  ro.r.text(mid, ro.r.par('usr')[2] * 2.5,
    srt=90, adj=1, labels=ro.r(labstr), xpd=ro.r('TRUE'), cex=cex_axis)

  # Add title
  if main:
    ro.r.title(main=ro.r(main), **{'cex.main' : cex_main})

  # Close PDF
  if outname:
    gr.dev_off()

  # Write table
  if taboutname:
    fh = open(taboutname, 'w')
    for optidx in range(len(fldopts)):
      if fldvals[optidx]:
        fh.write('%s\t%s\n' % (fldopts[optidx], fldvals[optidx]))
    fh.close()
  
def plotbool(boolinfo, rep, dv, cat, horiz=True, outname=None):
  
  if horiz:
    sortrev = False
    textoff = 0
    horopt = ro.r('TRUE')
    xlab = ro.r('"Proportion using procedure"')
    ylab = ''
    xlim = ro.FloatVector([0, 1])
    ylim = ro.r('NULL')
    legpos = 'bottomright'
  else:
    sortrev = True
    textoff = 2
    horopt = ro.r('FALSE')
    xlab = ''
    ylab = ro.r('"Proportion using procedure"')
    xlim = ro.r('NULL')
    ylim = ro.FloatVector([0, 1])
    legpos = 'topright'

  boolvals = [boolinfo[step][dv][cat] for step in boolsteps]
  boolzip = zip(boolsteps, boolvals)
  boolzip.sort(key=operator.itemgetter(1), reverse=sortrev)
  booldict = collections.OrderedDict(boolzip)
  
  xtk = arange(len(booldict))
  barcol = []
  bargrp = []
  for step in booldict:
    fstep = fullvar(step, rep[0].keys())
    fstep0 = fstep.split('-')[0]
    if fstep0 in grppal:
      barcol.append(grppal[fstep0])
      bargrp.append(grppal.keys().index(fstep0))
    else:
      barcol.append(grppal['def'])
      bargrp.append(len(grppal) - 1)
  
  rvec = ro.FloatVector(booldict.values())
  popts = {
    'xlab' : xlab,
    'ylab' : ylab,
    'xlim' : xlim,
    'ylim' : ylim,
    'cex.axis' : cex_axis,
    'cex.lab' : cex_lab,
    'horiz' : horopt,
    'col' : ro.StrVector(barcol)
  }

  # Open PDF
  if outname:
    gr.pdf(outname)

  # Set margins
  ro.r.par(mar=ro.FloatVector(mar))
  
  midpt = ro.r.barplot(rvec, **popts)
  a2z = [chr(x) for x in range(ord('a'), ord('z') + 1)]
  labs = booldict.keys()
  for labidx in range(len(labs)):
    if labs[labidx] == 'coplanar':
      labs[labidx] = 'aux'
  labstr = 'c(' + \
    ', '.join(['expression(phantom("%s")*"%s")' % (a2z, lab) for lab in labs]) + \
    ')'
  if horiz:
    ro.r.text(ro.r.par('usr')[textoff] - 0.025, midpt, 
      srt=0, adj=1, labels=ro.r(labstr), xpd=ro.r('TRUE'), cex=cex_axis)
  else:
    ro.r.text(midpt, ro.r.par('usr')[2] - 0.025, 
      srt=90, adj=1, labels=ro.r(labstr), xpd=ro.r('TRUE'), cex=cex_axis)
  
  # Add legend
  legvals = [grpname[grp] for grp in grpcol]
  ro.r.legend(legpos, ro.StrVector(legvals),
    fill=ro.StrVector(grppal.values()), cex=cex_leg, inset=0.025)
  #ro.r.legend('topright', ro.StrVector(legvals),
  #  fill=ro.StrVector(grppal.values()), cex=cex_leg)
  
  # Close PDF
  if outname:
    gr.dev_off()

def plotmissing(stepinfo, rep, usestrict, dv, steps, main=None, outname=None):
  
  missdict = collections.OrderedDict({})
  space = []
  col = []
  
  for grp in stepgrp:
    if steps[0] in stepgrp[grp]:
      break
  midcol = grpcol[grp]

  if midcol in divcols:
    colopts = ro.r('brewer.pal(%d, "%s")' % (len(steps) * 2, midcol))
    colopts = colopts[:len(steps)]
  else:
    colopts = ro.r('brewer.pal(%d, "%s")' % (len(steps) + 2, midcol))
    colopts = [colopts[x] for x in range(1, len(colopts) - 1)]

  for stepidx in range(len(steps)):
    
    step = steps[stepidx]

    stepvars = [v for v in rep[0].keys() if 
      re.search('^(\w+\-){0,}' + step + '\-', v) and
        not v.endswith('bool')]
    
    if 'strict' in stepinfo[step][stepvars[0]]:
      if usestrict:
        rule = 'strict'
      else:
        rule = 'lenient'
    else:
      rule = 'mand'
  
    missval = [1 - stepinfo[step][fld][rule][dv]['missing'] for fld in stepvars]
    misszip = zip(stepvars, missval)
    missdict.update(collections.OrderedDict(sorted(
      misszip, key=operator.itemgetter(1), reverse=True
    )))

    space.extend([0.8] + [0.2 for x in range(len(stepvars) - 1)])
    col.extend([colopts[stepidx] for x in range(len(stepvars))])

  barvec = ro.FloatVector(missdict.values())
  popts = {
    'ylim' : ro.FloatVector([0, 1]),
    'ylab' : ro.r('"Proportion describing parameter"'),
    'space' : ro.FloatVector(space),
    'col' : ro.StrVector(col),
    'cex.lab' : cex_lab,
    'cex.axis' : cex_axis
  }

  if outname:
    gr.pdf(outname, width=14)
  
  midpt = ro.r.barplot(barvec, **popts)
  
  a2z = [chr(x) for x in range(ord('a'), ord('z') + 1)]
  labstr = 'c(' + \
    ', '.join(['expression(phantom("%s")*"%s")' % (a2z, lab.split('-')[-1]) for lab in missdict]) + \
    ')'
  ro.r.text(midpt, ro.r.par('usr')[2] - 0.025, srt=90, adj=1, labels=ro.r(labstr), xpd=ro.r('TRUE'), cex=cex_axis)
  
  if main:
    ro.r.title(main=ro.r(main), **{'cex.main' : cex_main})
  
  if grp in legpos:
    pos = legpos[grp]
  else:
    pos = 'topright'
  legval = [stepname[step] for step in steps]
  ro.r.legend(pos, ro.StrVector(legval),
    fill=ro.StrVector(colopts), inset=0.025, cex=cex_leg)

  if outname:
    gr.dev_off()

def rulemand(dummy):
  return True
def makerulebool(step):
  return {
    'strict' : lambda r: r[step] == 'TRUE',
    'lenient' : lambda r: r[step] != 'FALSE'
}

def fullvar(var, vars):
  bstr = '^(\w+\-){0,}' + var + '\-bool$'
  matchvars = [v for v in vars if re.search(bstr, v)]
  if len(matchvars) == 0:
    return False
  elif len(matchvars) > 1:
    warnings.warn('too many matchvars: %s' % var)
  return matchvars[0]

def getincl(rep):
  incl = {}
  for step in mandsteps:
    incl[step] = { 'mand' : rulemand }
  for step in boolsteps:
    incl[step] = {}
    stepvar = fullvar(step, rep[0].keys())
    incl[step] = makerulebool(stepvar)
  return incl

def procbool(rep):
  
  repvars = rep[0].keys()
  incl = getincl(rep)
  boolinfo = {}
  
  for step in boolsteps:
    fstep = fullvar(step, repvars)
    vals = [r[fstep] for r in rep]
    uvals = list(set(vals))
    tmpcount = {}
    tmpprop = {}
    nincl = len([v for v in vals if v != 'n/a'])
    for uval in uvals:
      nval = float(len([v for v in vals if v == uval]))
      tmpcount[uval] = nval
      if uval != 'n/a':
        tmpprop[uval] = nval / nincl
    
    boolinfo[step] = {
      'count' : collections.OrderedDict(
        sorted(tmpcount.iteritems(), key=operator.itemgetter(1),
        reverse=True
      )),
      'prop' : collections.OrderedDict(
        sorted(tmpprop.iteritems(), key=operator.itemgetter(1),
        reverse=True
      ))
    }

  return boolinfo

def procsteps(rep):
  
  repvars = rep[0].keys()
  incl = getincl(rep)
  
  stepct = {}

  for step in allsteps:
    
    stepct[step] = {}
    stepvars = [v for v in repvars if 
      re.search('^(\w+\-){0,}' + step + '\-', v) and
      not v.endswith('bool')]
    
    for var in stepvars:
      
      stepct[step][var] = {}
      
      for rule in incl[step]:
      
        increc = [r for r in rep if incl[step][rule](r)]
        stepct[step][var][rule] = {}
        
        vals = [str(r[var]).lower().split('; ') for r in increc]
        #vals = [r[var].lower().split('; ') for r in increc]
        fvals = reduce(operator.add, vals)
        uvals = list(set(fvals))
        
        tmpcount = {'missing' : 0}
        tmpprop = {'missing' : 0}
        tmpgprop = {}
        
        nincl1 = float(len([v for v in vals if 'n/a' not in v]))
        nincl2 = float(len([v for v in vals if 
          'n/a' not in v and 'missing' not in v]))
        #nincl1 = float(len([v for v in fvals if v != 'n/a']))
        #nincl2 = float(len([v for v in fvals if 
        #  v not in ['n/a', 'missing']]))
        
        for uval in uvals:
          
          #nval = len([v for v in fvals if v == uval])
          nval = len([v for v in vals if uval in v])
          tmpcount[uval] = nval
          if uval != 'n/a':
            tmpprop[uval] = nval / nincl1
            if uval != 'missing':
              tmpgprop[uval] = nval / nincl2
        
        #stepct[step][var][rule] = {
        #  'count' : collections.OrderedDict(sorted(
        #    tmpcount.iteritems(), key=operator.itemgetter(1), 
        #    reverse=True
        #  )),
        #  'prop' : collections.OrderedDict(sorted(
        #    tmpprop.iteritems(), key=operator.itemgetter(1),
        #    reverse=True
        #  )),
        #  'gprop' : collections.OrderedDict(sorted(
        #    tmpgprop.iteritems(), key=operator.itemgetter(1),
        #    reverse=True
        #  ))
        #}
        stepct[step][var][rule] = {
          'count' : collections.OrderedDict(sorted(
            tmpcount.iteritems(), cmp=lambda x, y: cmp( (-x[1], x[0]), (-y[1], y[0]) )
          )),
          'prop' : collections.OrderedDict(sorted(
            tmpprop.iteritems(), cmp=lambda x, y: cmp( (-x[1], x[0]), (-y[1], y[0]) )
          )),
          'gprop' : collections.OrderedDict(sorted(
            tmpgprop.iteritems(), cmp=lambda x, y: cmp( (-x[1], x[0]), (-y[1], y[0]) )
          ))
        }

  return stepct

def getmcc(rep):
  
  missmc = [r for r in rep if r['mccorrectbool'] == 'missing']
  truemc = [r for r in rep if r['mccorrectbool'] == 'TRUE']
  falsemc = [r for r in rep if r['mccorrectbool'] == 'FALSE']
  
  print len(truemc), len(falsemc), len(missmc)

def plotext(rep, exttype='vol', outname=None):
  
  if exttype == 'vox':
    ext = [r['mod-gmod-mccorrect-ext'] for r in rep]
    ext = [x for x in ext if x not in ['missing', 'n/a']
      and x.endswith('voxels')]
    ext = [int(re.sub('\s+voxels', '', x)) for x in ext]
    cutoff = 250
    skip = 2
    xlab='Cluster extent threshold (# voxels)',
  elif exttype == 'vol':
    ext = [r['misc-extvol'] for r in rep]
    ext = [x for x in ext if x not in ['missing', 'n/a']]
    ext = [float(x) for x in ext]
    cutoff = 3000
    skip = 2
    xlab=ro.r('expression(paste("Cluster extent threshold (", mm^3, ")"))')
  
  print sorted(ext)
  print 'Median cluster extent: %f' % (median(ext))
  print 'Modal cluster extent: %f' % (scipy.stats.stats.mode(ext)[0][0])
  
  #histgap(
  #  {'ext' : ext},
  #  cutoff=cutoff,
  #  xlab='Cluster extent threshold (# voxels)',
  #  ylab='Count',
  #  skip=skip,
  #  main='expression(bold("Inference: Extent"))',
  #  outname=outname
  #)
  histgap(
    {'ext' : ext},
    cutoff=cutoff,
    nbins=8,
    xlab=xlab,
    ylab='Count',
    skip=skip,
    main='expression(bold("Inference: Extent"))',
    outname=outname
  )

def plotsmooth(rep, outname=None):
  
  fwhm = [float(r['proc-smooth-kernel'].split('mm')[0]) for r in rep
    if r['proc-smooth-bool'] == 'TRUE' and
    r['proc-smooth-kernel'].find('fwhm') != -1]
  breaks = ro.FloatVector(arange(floor(min(fwhm)) - 0.5, ceil(max(fwhm)) + 1))
  
  histgap(
    {'fwhm' : fwhm},
    cutoff=None,
    breaks=breaks,
    xlab='Full width at half maximum (mm)',
    ylab='Count',
    main='expression(bold("Spatial smoothing"))',
    outname=outname
  )

def procsoft(rep):
  
  softpck = []
  softver = []
  softcom = []

  for r in rep:
    pck = r['softpck'].split('; ')
    ver = r['softver'].split('; ')
    com = ['%s %s' % (pck[sidx], ver[sidx]) 
      for sidx in range(len(pck))]
    softpck.append(pck)
    softver.append(ver)
    softcom.append(com)
  
  flatpck = reduce(operator.add, softpck)
  print set(flatpck)
  return softcom

def getpwr(rep, esize, sig):
    
  alt = 'greater'
  #alt = 'two.sided'

  for r in rep:
    
    if 'subjinfo' not in r:
      continue
    
    if r['subjinfo']['ngrp'] == 1:

      n = r['subjinfo']['subjtot']
      if n == 1:
        continue
      
      pwrtmp = ro.r('pwr.t.test(n=%d, d=%f, sig.level=%f, alternative="%s")' % 
        (n, esize, sig, alt)
      )
      r['subjinfo']['pwr'] = pwrtmp[3][0]

    elif r['subjinfo']['ngrp'] == 2:
      
      ns = [r['subjinfo']['grp'][g] for g in r['subjinfo']['grp']]
      if 1 in ns:
        continue
      pwrtmp = ro.r('pwr.t2n.test(n1=%d, n2=%d, d=%f, sig.level=%f, alternative="%s")' %
        (ns[0], ns[1], esize, sig, alt)
      )
      r['subjinfo']['pwr'] = pwrtmp[4][0]

def getcumpwr(rep, esize, sig):
  
  subjinfo(rep)
  getpwr(rep, esize, sig)
  
  pwrrep = [r for r in rep if 'subjinfo' in r
    and r['subjinfo']['gt1']]

  pwr1 = [100 * r['subjinfo']['pwr'] for r in pwrrep
    if r['subjinfo']['ngrp'] == 1]
  pwr2 = [100 * r['subjinfo']['pwr'] for r in pwrrep
    if r['subjinfo']['ngrp'] == 2]

  cumpwr1 = []
  cumpwr2 = []

  pwrstep = 0.1
  pwrbins = arange(0, 100 + pwrstep, pwrstep)

  for pwrval in pwrbins:
    nstu1 = len([r for r in pwr1 if r >= pwrval])
    nstu2 = len([r for r in pwr2 if r >= pwrval])
    cumpwr1.append(nstu1 / float(len(pwr1)))
    cumpwr2.append(nstu2 / float(len(pwr2)))
  
  nrep1 = len([r for r in pwr1 if r >= 80])
  nrep2 = len([r for r in pwr2 if r >= 80])
  print 'Effect size %f, threshold %f' % (esize, sig)
  print 'One-group: %f' % (nrep1 / float(len(pwr1))), mean(pwr1)
  print 'Two-group: %f' % (nrep2 / float(len(pwr2))), mean(pwr2)

  return pwrbins, cumpwr1, cumpwr2

def plotpwr(rep, esize, sig, nbins=10, binsize=None, outname=None):

  getpwr(rep, esize, sig)
  
  pwrrep = [r for r in rep if 'subjinfo' in r 
    and r['subjinfo']['gt1']]

  pwr1 = [r['subjinfo']['pwr'] for r in pwrrep
    if r['subjinfo']['ngrp'] == 1]
  pwr2 = [r['subjinfo']['pwr'] for r in pwrrep
    if r['subjinfo']['ngrp'] == 2]

  print float(len([p for p in pwr1 if p >= 0.8])) / len(pwr1)
  print float(len([p for p in pwr2 if p >= 0.8])) / len(pwr2)
  
  if binsize:
    pwrbin = arange(0, 1 + binsize, binsize)
  else:
    pwrhist = ro.r.hist(ro.FloatVector(pwr1 + pwr2), nbins, plot=False)
    pwrbin = list(pwrhist[0])
  binsize = pwrbin[1] - pwrbin[0]

  histgap(
    collections.OrderedDict((
      ('One group', pwr1),
      ('Two groups', pwr2)
    )),
    maxbin=binsize + 1,
    nbins=nbins,
    binsize=binsize,
    breaks=pwrbin,
    cutoff=None,
    main='expression(bold(paste("Power at ", bolditalic("d"), "=%0.1f, ", bolditalic("p"), "=%0.3f")))' % (esize, sig),
    xlab='Power (%)',
    ylab='Count',
    skip=2,
    outname=outname
  )
  return
  if binsize:
    pwrbin = arange(0, 1 + binsize, binsize)
  else:
    pwrhist = ro.r.hist(ro.FloatVector(pwr1 + pwr2), nbins, plot=False)
    pwrbin = list(pwrhist[0])
  binsize = pwrbin[1] - pwrbin[0]

  pwr1ct = list(ro.r.hist(ro.FloatVector(pwr1), 
    breaks=ro.FloatVector(pwrbin), plot=False)[1])
  pwr2ct = list(ro.r.hist(ro.FloatVector(pwr2), 
    breaks=ro.FloatVector(pwrbin), plot=False)[1])

  # Open PDF
  if outname:
    gr.pdf(outname)
  
  pwr1str = ', '.join([str(p) for p in pwr1ct])
  pwr2str = ', '.join([str(p) for p in pwr2ct])
  borstr = ', '.join(
    ['"%s"' % 
      {True : 'black', False : 'white'}[pwr1ct[pidx] > 0 or pwr2ct[pidx] > 0] 
      for pidx in range(len(pwr1ct))
    ]
  )
  s = """modbarplot(rbind(c(%s), c(%s)), 
    beside=FALSE,
    border=c(%s),
    main=expression(paste("Power at ", italic("d"), "=%0.1f, ", italic("p"), "=%0.3f")),
    xlab="Power (%%)",
    ylab="Count")""" \
    % (pwr1str, pwr2str, borstr, esize, sig)
  
  mid = ro.r(s)
  
  lab = ['%0.1f' % ((pwrbin[aidx] + binsize / 2) * 100) for aidx in range(len(pwrbin) - 1)]
  
  # Draw axis
  axopts = {
    'labels' : ro.StrVector(lab),
    'las' : 2
  }
  ro.r.axis(1, at=mid, **axopts)

  # Add legend
  if pwr1ct[0] > pwr1ct[-1]:
    legpos = 'topright'
  else:
    legpos = 'topleft'
  ro.r.legend(
    legpos,
    ro.StrVector(('One Group', 'Two Groups')),
    fill=ro.r('gray.colors(2)'),
    inset=0.025
  )
  
  # Close PDF
  if outname:
    gr.dev_off()

def plothpf(rep, cutoff, nbins=10, binsize=None, outname=None):
  
  filtval = [r['mod-smod-filter-filtband'] for r in rep if 
    r['mod-smod-filter-filttype'] == 'high-pass'
    and r['des-edes-destype'] == 'blocked'
  ]
  
  for r in rep:
    
    if r['mod-smod-filter-filttype'] != 'high-pass':
      continue

    val = r['mod-smod-filter-filtband']

    secsrc = re.search('^(\d+)s$', val)
    if secsrc:
      r['hpfsec'] = float(secsrc.groups()[0])
    
    hzsrc = re.search('hz$', val, re.I)
    if hzsrc:
      hzstrip = re.sub('\s*hz', '', val, flags=re.I)
      r['hpfsec'] = 1 / float(hzstrip)

  evthpf = [r['hpfsec'] for r in rep if 
    'hpfsec' in r and r['des-edes-destype'] == 'event-related']
  blkhpf = [r['hpfsec'] for r in rep if 
    'hpfsec' in r and r['des-edes-destype'] == 'blocked']
  othhpf = [r['hpfsec'] for r in rep if 
    'hpfsec' in r and r['des-edes-destype'] in ['n/a', 'missing']]

  print median(evthpf)
  print median(blkhpf)
  print median(othhpf)

  print median(evthpf + blkhpf + othhpf)

  print(sorted(evthpf + blkhpf + othhpf))
  
  histgap(
    collections.OrderedDict((
      ('Event', evthpf),
      ('Block', blkhpf),
      ('Other', othhpf)
    )),
    cutoff=cutoff,
    nbins=nbins,
    binsize=binsize,
    xlab='Filter cutoff (s)',
    ylab='Count',
    main='expression(bold("Temporal filtering"))',
    skip=2,
    outname=outname
  )

def histgap(histdat, cutoff, nbins=10, binsize=None, maxbin=None, breaks=None, main='', xlab='', ylab='', skip=1, outname=None):
  
  grps = histdat.keys()
  trimdat = {}
  xtradat = {}
  
  if not cutoff:
    cutoff = inf

  for grp in grps:

    trimdat[grp] = [dat for dat in histdat[grp] if dat < cutoff]
    xtradat[grp] = [dat for dat in histdat[grp] if dat >= cutoff]
  
  compall = reduce(operator.add, histdat.values())
  trimall = reduce(operator.add, trimdat.values())
  xtraall = reduce(operator.add, xtradat.values())
  
  if not maxbin:
    maxbin = max(trimall)
  
  if binsize:
    gapbin = arange(0, maxbin, binsize)
  else:
    gaphist = ro.r.hist(ro.FloatVector(trimall), nbins, plot=False)
    gapbin = list(gaphist[0])
  binsize = gapbin[1] - gapbin[0]
  binxtra = arange(float(gapbin[0]), float(max(compall) + binsize), float(binsize))

  if breaks == None:
    breaks = binxtra
  
  histct = {}
  trimct = {}

  for grp in grps:

    histct[grp] = list(
      ro.r.hist(
        ro.FloatVector(histdat[grp]),
        breaks=ro.FloatVector(breaks),
        plot=False
      )[1]
    )
    print len(gapbin), len(histct[grp])
    trimct[grp] = [histct[grp][ctidx] for ctidx in range(len(histct[grp]))]
    #trimct[grp] = [histct[grp][ctidx] for ctidx in range(len(gapbin))]
    #    breaks=ro.FloatVector(binxtra),
  
  ctxtra = ro.r.hist(
    ro.FloatVector(compall), 
    breaks=ro.FloatVector(breaks),
    plot=False
  )[1]
  #  breaks=ro.FloatVector(binxtra),
  
  if not isinf(cutoff):

    gappos = int(floor(float(cutoff) / binsize))
    incidx = []
    
    for ctidx in range(gappos):
      
      break2cut = ctxtra[ctidx:gappos]
      if sum(break2cut) > 0:
        incidx.append(ctidx)
      else:
        break

    gapidx = ctidx + 2
    incidx.extend([ctidx] * 2)
    
    firstpast = False
    for ctidx in range(gappos, len(ctxtra)):
      cut2break = ctxtra[gappos:ctidx+1]
      if sum(cut2break) > 0:
        if not firstpast:
          firstpast = True
          incidx.append(ctidx - 1)
        incidx.append(ctidx)

    incidx.append(ctidx + 1)
    
    ctxtra = list(ctxtra)
    ctxtra.append(0)
    for grp in grps:
      histct[grp].append(0)

  else:

    incidx = range(len(ctxtra))
    gapidx = inf

  borcol = [{True : '"black"', False : '"white"'}[ctxtra[idx] > 0] for idx in incidx]
  borstr = ', '.join(borcol)
  
  cstr = []
  for grp in grps:
    cstr.append(
      ', '.join([str(histct[grp][idx]) for idx in incidx])
    )
  cstr.reverse()
  comstr = ', '.join(['c(%s)' % s for s in cstr])
  if len(histdat) > 1:
    datstr = 'rbind(%s)' % (comstr)
  else:
    datstr = comstr
  
  s0 = """modbarplot(
    %s,
    border=c(%s), 
    cex.lab=%f,
    cex.axis=%f,
    beside=FALSE)""" \
    % (datstr, borstr, cex_lab, cex_axis)

  ro.r(s0)
  
  ticks = ro.r.axTicks(2)

  # Open PDF
  if outname:
    gr.pdf(outname)
  
  #ro.r.par(mar=mar)
  # Set margins
  ro.r.par(mar=ro.FloatVector(mar))

  ymax = max(ctxtra)
  if ticks[-1] < ymax:
    print ticks
    ymax = ticks[-1] + (ticks[-1] - ticks[-2])
    ylim = 'c(%f, %f)' % (ymax * -0.01, ymax)
  else:
    ylim = 'NULL'

  s = """modbarplot(
    %s,
    border=c(%s), 
    ylim=%s,
    cex.lab=%f,
    cex.axis=%f,
    beside=FALSE)""" \
    % (datstr, borstr, ylim, cex_lab, cex_axis)
  mid = ro.r(s)
  
  #
  at = [mid[aidx] for aidx in range(len(mid)) if aidx != gapidx - 1]
  lab = [breaks[incidx[aidx]] + binsize / 2 for aidx in range(len(incidx)) if aidx != gapidx - 1]
  
  # 
  #at = [at[idx] for idx in range(len(at)) if idx % skip == 0]
  #lab = [lab[idx] for idx in range(len(lab)) if idx % skip == 0]
  for labidx in range(len(lab)):
    if labidx + 1 < gapidx:
      modidx = labidx
    else:
      modidx = labidx - gapidx
    if modidx % skip != 0:
      lab[labidx] = ''
    elif lab[labidx].is_integer():
      lab[labidx] = '%d' % (lab[labidx])
    else:
      lab[labidx] = '%0.1f' % (lab[labidx] * 100)
  print 'lab', lab
  print 'at', at

  # Set label offsets
  ro.r('par(mgp=c(%f, %f, %f))' % (mgpx[0], mgpx[1], mgpx[2]))
  ro.r.title(xlab=xlab, **{'cex.lab' : cex_lab})
  ro.r('par(mgp=c(%f, %f, %f))' % (mgpy[0], mgpy[1], mgpy[2]))
  ro.r.title(ylab=ylab, **{'cex.lab' : cex_lab})
  if main:
    ro.r.title(main=ro.r(main), **{'cex.main' : cex_main})

  # Draw axis
  axopts = {
    'labels' : ro.StrVector(lab),
    'cex.axis' : cex_axis,
    'las' : 2
  }
  ro.r.axis(1, at=ro.FloatVector(at), **axopts)
  
  # Add axis break
  if not isinf(cutoff):
    ro.r('axis.break(1, breakpos=%f, brw=0.03)' % (mid[gapidx - 1]))
  
  # 
  if ctxtra[0] > ctxtra[-1]:
    legpos = 'topright'
  else:
    legpos = 'topleft'
  
  # Add legend
  if len(histdat) > 1:
    ro.r.legend(
      legpos,
      ro.StrVector(ro.StrVector(histdat.keys())),
      fill=ro.r('gray.colors(%d, start=0.9, end=0.3)' % (len(histdat))),
      cex=cex_leg,
      inset=0.025
    )

  # Close PDF
  if outname:
    gr.dev_off()
  
  
def plotsubjs(rep, cutoff, nbins=10, binsize=None, outname=None):
  
  pwrrep = [r for r in rep if 'subjinfo' in r 
    and r['subjinfo']['gt1']]

  ss1 = [r['subjinfo']['subjavg'] for r in pwrrep
    if r['subjinfo']['ngrp'] == 1]
  ss2 = [r['subjinfo']['subjavg'] for r in pwrrep
    if r['subjinfo']['ngrp'] == 2]
  
  print len(pwrrep), len(ss1), len(ss2)
  print min(ss1), max(ss1), median(ss1)
  print min(ss2), max(ss2), median(ss2)

  histgap(
    collections.OrderedDict((
      ('One group', ss1),
      ('Two groups', ss2)
    )),
    cutoff=cutoff,
    nbins=nbins,
    binsize=binsize,
    main='expression(bold("Sample size"))',
    xlab='Subjects per group',
    ylab='Count',
    skip=2,
    outname=outname)
  return
  
  
  ss1trim = [ss for ss in ss1 if ss < cutoff]
  ss2trim = [ss for ss in ss2 if ss < cutoff]
  ss1xtra = [ss for ss in ss1 if ss >= cutoff]
  ss2xtra = [ss for ss in ss2 if ss >= cutoff]
  
  if binsize:
    ssbin = range(0, int(max(ss1trim + ss2trim)), binsize)
  else:
    sshist = ro.r.hist(ro.FloatVector(ss1trim + ss2trim), nbins, plot=False)
    ssbin = list(sshist[0])
  binsize = ssbin[1] - ssbin[0]
  ssbinxtra = range(int(ssbin[0]), int(max(ss1 + ss2) + binsize), int(binsize))
  maxbintrim = ssbin[-1]
  
  ss1ct = list(ro.r.hist(ro.FloatVector(ss1), breaks=ro.FloatVector(ssbinxtra), plot=False)[1])
  ss1cttrim = [ss1ct[ssidx] for ssidx in range(len(ssbin))]
  
  ss2ct = list(ro.r.hist(ro.FloatVector(ss2), breaks=ro.FloatVector(ssbinxtra), plot=False)[1])
  ss2cttrim = [ss2ct[ssidx] for ssidx in range(len(ssbin))]
  
  ssctxtra = ro.r.hist(
    ro.FloatVector(ss1 + ss2), 
    breaks=ro.FloatVector(ssbinxtra),
    plot=False
  )[1]
  minbinxtra = min([ssbinxtra[ssidx] for ssidx in range(len(ssbinxtra))
    if ssidx > len(ssbin) and ssidx < len(ssctxtra) and ssctxtra[ssidx] > 0])
  
  if len(ss1xtra + ss2xtra) > 0:
    ssbreak = len(ssbin)
    ssbin.extend([0] * 1 + [minbinxtra - binsize])
    ss1cttrim.extend([0] * 2)
    ss2cttrim.extend([0] * 2)

  for ssidx in range(len(ssbin), len(ssbinxtra)):
    
    if ssbinxtra[ssidx] < ssbin[-1]:
      continue

    if ssidx >= len(ssctxtra):
      continue
    
    if ssctxtra[ssidx] == 0:
      continue
    
    ssbin.append(ssbinxtra[ssidx])
    if ssidx < len(ss1ct):
      ss1cttrim.append(ss1ct[ssidx])
    else:
      ss1cttrim.append(0)
    if ssidx < len(ss2ct):
      ss2cttrim.append(ss2ct[ssidx])
    else:
      ss2cttrim.append(0)
    
  borcol = []
  for ssidx in range(len(ssbin)):
    if ss1cttrim[ssidx] > 0 or ss2cttrim[ssidx] > 0:
      borcol.append('black')
    else:
      borcol.append('white')
  
  ssbin.append(ssbin[-1] + binsize)
  ss1cttrim.append(0)
  ss2cttrim.append(0)
  borcol.append('white')
  borcol = ro.StrVector(borcol)
  baropts = {
    'cex.names' : 0.75,
    'beside' : ro.r('FALSE'),
    'border' : borcol
  }
  
  # Open PDF
  if outname:
    gr.pdf(outname)
  
  ss1str = ', '.join([str(ss) for ss in ss1cttrim])
  ss2str = ', '.join([str(ss) for ss in ss2cttrim])
  borstr = ', '.join(['"%s"' % b for b in borcol])
  s = """rbind(c(%s), c(%s)), beside=TRUE, border=c(%s)""" \
    % (ss1str, ss2str, borstr)
  mid = ro.r("""modbarplot(rbind(c(%s), c(%s)), 
    beside=FALSE,
    xlab="Subjects per Group",
    ylab="Count",
    border=c(%s))""" \
    % (ss1str, ss2str, borstr))
  
  at = [mid[aidx] for aidx in range(len(mid)) if aidx != ssbreak]
  lab = [ssbin[aidx] + binsize / 2 for aidx in range(len(ssbin)) if aidx != ssbreak]
  
  # Draw axis
  axopts = {
    'labels' : ro.FloatVector(lab),
    'las' : 2
  }
  ro.r.axis(1, at=ro.FloatVector(at), **axopts)
  
  # Add axis break
  ro.r('axis.break(1, breakpos=%f, brw=0.03)' % (mid[ssbreak]))

  # Add legend
  ro.r.legend(
    'topright', 
    ro.StrVector(('One Group', 'Two Groups')),
    fill=ro.r('gray.colors(2)'),
    inset=0.025
  )
  
  # Close PDF
  if outname:
    gr.dev_off()

def subjinfo(rep):
  
  for r in rep:
    if r['des-subjs-nsubj'] == 'missing':
      continue
    tmpdict = {'grp' : {}}
    tmpdict['gt1'] = True
    if r['des-subjs-nsubj'].find(';') == -1:
      tmpdict['grp'] = {'all' : int(r['des-subjs-nsubj'])}
      if int(r['des-subjs-nsubj']) == 1:
        tmpdict['gt1'] = False
      if r['des-subjs-nfemale'] != 'missing':
        tmpdict['nfem'] = {'all' : int(r['des-subjs-nfemale'])}
    else:
      r['xnsubj'] = {}
      for grp in r['des-subjs-nsubj'].split(';'):
        sgrp = grp.split(':')
        grpname = sgrp[0].split('group')[0].strip()
        tmpdict['grp'][grpname] = int(sgrp[1])
        if int(sgrp[1]) == 1:
          tmpdict['gt1'] = False
      if r['des-subjs-nfemale'] != 'missing':
        tmpdict['nfem'] = {}
        for grp in r['des-subjs-nfemale'].split(';'):
          sgrp = grp.split(':')
          grpname = sgrp[0].split('group')[0].strip()
          tmpdict['nfem'][grpname] = int(sgrp[1])

    tmpdict['ngrp'] = len(tmpdict['grp'])

    tmpdict['subjtot'] = sum([tmpdict['grp'][g] for g in tmpdict['grp']])
    tmpdict['subjavg'] = mean([tmpdict['grp'][g] for g in tmpdict['grp']])

    if r['des-subjs-nfemale'] != 'missing':
      tmpdict['propfem'] = mean([ 
        tmpdict['nfem'][g] / float(tmpdict['grp'][g]) 
        for g in tmpdict['grp']
        ])

    r['subjinfo'] = tmpdict
  
  return
  nsubjws = [r['subjinfo']['subjavg'] for r in rep \
    if r['nsubj'] != 'missing' and r['subjinfo']['ngrp'] == 1 \
    and r['casestudy'] == 'FALSE']
  nsubjbs = [r['subjinfo']['subjavg'] for r in rep \
    if r['nsubj'] != 'missing' and r['subjinfo']['ngrp'] > 1 \
    and r['casestudy'] == 'FALSE']
  print median(nsubjws), median(nsubjbs)

  pfemws = [r['subjinfo']['propfem'] for r in rep \
    if r['nsubj'] != 'missing' and r['nfemale'] != 'missing' \
    and r['subjinfo']['ngrp'] == 1 and r['casestudy'] == 'FALSE' ]
  pfembs = [r['subjinfo']['propfem'] for r in rep \
    if r['nsubj'] != 'missing' and r['nfemale'] != 'missing' \
    and r['subjinfo']['ngrp'] > 1 and r['casestudy'] == 'FALSE' ]
  print mean(pfemws), mean(pfembs)

def stepord(rep, steps):
  
  ordct = {}

  for r in rep:
    
    steplist = r['steplist'].split(', ')
    stepidx = [(step, steplist.index(step)) for step in steps
      if step in steplist]
    if len(stepidx) != len(steps):
      continue
    stepidx.sort(key=operator.itemgetter(1))
    stepidx = [step[0] for step in stepidx]
    
    stepstr = '-'.join(stepidx)

    if stepstr not in ordct:
      ordct[stepstr] = 0
    ordct[stepstr] += 1

  return ordct

def addboolstream(rep, skipdispl=True):
  
  ustream = []

  for r in rep:
    
    bstmp = collections.OrderedDict({})

    for step in boolsteps:

      if skipdispl and step in stepgrp['displ']:
        continue

      fstep = fullvar(step, rep[0].keys())
      if r[fstep] == 'TRUE':
        bstmp[step] = True
      else:
        bstmp[step] = False

    r['boolstream'] = bstmp

    if bstmp not in ustream:
      ustream.append(bstmp)

  print len(ustream), len(rep)
  
def mergerep(rep):
  
  for r in rep:
    for fld in r:
      val = r[fld]
      if type(val) == str:
        val = r[fld].lower()
      if fld in merge and val in merge[fld]:
        r[fld] = merge[fld][val]

  return rep

def augrep(rep):
  
  steplists = [r['steplist'].split(', ') for r in rep]
  stepcat = reduce(operator.add, steplists)
  usteps = list(set(stepcat))
  usteps = [step for step in usteps if 
    not fullvar(step, rep[0].keys())]

  for ridx in range(len(rep)):
    
    for step in usteps:
      stepname = 'proc-%s-bool' % (step)
      if step in steplists[ridx]:
        rep[ridx][stepname] = 'TRUE'
      else:
        rep[ridx][stepname] = 'missing'

    if 'movement params' in rep[ridx]['mod-smod-regress'].split('; '):
      rep[ridx]['mod-smod-motreg-bool'] = 'TRUE'
    else:
      rep[ridx]['mod-smod-motreg-bool'] = 'missing'
    
    jtitle = rep[ridx]['jtitle'].lower()
    jtitle = re.sub('[\.;:].*', '', jtitle)
    jtitle = re.sub('\(.*', '', jtitle)
    rep[ridx]['misc-jtitle'] = jtitle.strip()

    softpck = rep[ridx]['misc-softpck'].split('; ')
    softver = rep[ridx]['misc-softver'].split('; ')
    softcom = []
    
    for softidx in range(len(softpck)):
      pck = softpck[softidx]
      ver = softver[softidx]
      if ver == 'missing':
        ver = '?'
      if pck == 'missing':
        com = 'missing'
      elif pck == 'custom':
        com = 'custom'
      else:
        com = '%s %s' % (softpck[softidx], ver)
      softcom.append(com)
    softstr = '; '.join(softcom)
    rep[ridx]['misc-softcom'] = softstr
    
    # Get extent threshold volume
    rep[ridx]['misc-extvol'] = getextvol(rep[ridx])

    # 
    if rep[ridx]['mod-smod-filter-bool'] != 'TRUE':
      rep[ridx]['mod-smod-filter-filttype'] = 'n/a'
      rep[ridx]['mod-smod-filter-filtband'] = 'n/a'

    if rep[ridx]['mod-smod-acf-bool'] != 'TRUE':
      rep[ridx]['mod-smod-acf-acfmethod'] = 'n/a'

    del rep[ridx]['mod-smod-type']

  return rep

def getextvol(rep):
  
  ext = rep['mod-gmod-mccorrect-ext']
  
  if ext == 'n/a':
    return 'n/a'

  if ext in ['missing', 'misisng']:
    return 'missing'

  ccmatch = re.search('(\d+)\s(mm3|ul)', ext, re.I)
  if ccmatch:
    try:
      return int(ccmatch.groups()[0])
    except:
      pass
  
  evmatch = re.search('(\d+)\svoxels$', ext)
  if not evmatch:
    return 'missing'
  try:
    nvox = int(evmatch.groups()[0])
  except:
    return 'missing'

  finalres = rep['misc-finalres']
  if not finalres or finalres in ['missing', 'misisng']:
    return 'missing'
  resdims = finalres.split('x')
  if len(resdims) != 3:
    return 'missing'
  try:
    voxsize = prod([float(dim) for dim in resdims])
  except:
    return 'missing'

  return nvox * voxsize

def readrep(repfile):
  
  replines = open(repfile, 'r').readlines()
  replines = [line.strip().split('\t') for line in replines]

  repcatvars = replines[0]
  repvars = replines[1]
  repdata = replines[2:]
  
  repdict = []
  for repline in repdata:
    tmpdict = dict([ \
      (repvars[idx], repline[idx]) \
      for idx in range(len(repline)) \
      if repvars[idx] != '' \
      ])
    if 'include' in tmpdict and tmpdict['include'] == 'yes':
      repdict.append(tmpdict)
  
  return repdict
