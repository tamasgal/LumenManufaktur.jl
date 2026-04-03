# API

```@index
```

## Parameters

```@docs
LMParameters
PDFIntegrationPoints
PMTModel
```

## PDF Tables

```@docs
PDFTableConfig
DirectMuonPDFTable
makepdf
save_pdftable
load_pdftable
MuonPDFTableConfig
MuonPDFTable
makemuonpdf
save_muonpdftable
load_muonpdftable
EMShowerPDFTableConfig
EMShowerPDFTable
makeemshowerpdf
save_emshowerpdftable
load_emshowerpdftable
BrightPointPDFTableConfig
BrightPointPDFTable
makebrightpointpdf
save_brightpointpdftable
load_brightpointpdftable
```

## Muon light

```@docs
directlightfrommuon
scatteredlightfrommuon
```

## Delta-ray light

```@docs
directlightfromdeltarays
scatteredlightfromdeltarays
deltarayenergyloss
```

## EM shower light

```@docs
directlightfromEMshower
scatteredlightfromEMshower
directlightfromEMshowers
scatteredlightfromEMshowers
```

## Bright-point light

```@docs
directlightfrombrightpoint
scatteredlightfrombrightpoint
```

## Combined light

```@docs
lightfrommuon
lightfromEMshower
lightfrombrightpoint
```

## Dispersion

```@docs
DispersionModel
BaileyDispersion
DispersionORCA
DispersionARCA
refractionindexphase
refractionindexgroup
```

## Absorption

```@docs
AbsorptionModel
DefaultAbsorption
absorptionlength
```

## Scattering

```@docs
ScatteringModel
Kopelevich
scatteringlength
ScatteringProbabilityModel
Scatteringp00075
scatteringprobability
henyey_greenstein
rayleigh
```

## EM shower profile

```@docs
Geanz
geanz
getprobability
getintegral
getlength
getmaximum
geant
```
