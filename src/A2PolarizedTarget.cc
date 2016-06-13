#include "A2PolarizedTarget.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Sphere.hh"
#include "G4Cons.hh"
#include "G4SubtractionSolid.hh"
#include "G4PVPlacement.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4VisAttributes.hh"
#include "A2MagneticField.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4Paraboloid.hh"
#include "G4Polycone.hh"

// Sam Abernethy, June 2016
// Changing A2PolarizedTarget to include Active Target
// Old version of A2PolarizedTarget.cc is in 'A2NotActivePolarizedTarget.cc'
// Implementing boolean 'active' as true or false based on input target material
// setting up (active == true) as standalone section
// setting up (active == false) with kapton cell, butanol cell, Cu cylinders but not yet helium between cylinders!
// helium between cylinders will be important

/* ********************************* STRUCTURE OF A2POLARIZEDTARGET.CC ******************************
 * Magnetic field is defined in fMagneticField, and SetMagneticField void function.
 * Mother Volume constructed.
 *
 * All other Volumes constructed:
 * Define length, radius, thickness (or x, y, z).
 * Define G4Tubs/G4Box.
 * Define G4LogicalVolume.
 * Define G4PVPlacement with use of G4ThreeVector.
 * SetVisAttributes.
 *
 * TypeMagneticCoils from DetectorSetup.mac can either be Solenoidal or Saddle:
 * If Solenoidal, construct cylinders.
 * If Saddle, construct boxes (and subtract them).
 *
 * Target cell, butanol, is defined (and hardcoded).
 * */
/* ************************************** VOLUMES AND MATERIALS *************************************
 * "Mother" Air volume -- fMyLogic
 * Outer SS cylinder -- SSOLogic
 * Outer Cu cylinder -- CUOLogic
 * Solenoidal magnetic coils --  NbTiCLogic, CUCLogic, EPCLogic
 * Saddle magnetic coils -- logicSaddleCoilsLayer1, logicSaddleCoilsLayer2
 * Middle Cu cylinder -- CUMALogic, CUMBLogic, CUMCLogic
 * Inner SS cylinder -- SSIALogic, SSIBLogic
 * Inner Cu cylinder -- CUIALogic, CUIBLogic, CUICLogic, CUIDLogic
 * Kapton cell -- KAPALogic, KAPBLogic, KAPCLogic, KAPDLogic
 * Butanol Target cell -- BTRGTLogic
 * He between inner SS and Cu cylinders -- HEALogic, ..., HEFLogic
 * Outer Ti window -- TIOWLogic
 * Cu ring -- CUOWLogic
 * Outer Al window -- ALOWLogic
 * Inner Al window -- ALIWLogic
 * Cu bit between coils and inner Al window -- CUBALogic, CUBBLogic
 * Middle Ti window -- TIMWLogic
 * SS ring -- SSIWLogic
 * Inner Ti window -- TIIWLogic
 *
 * l = length, r = radius, t = thickness
 * distances are HALF-LENGTHS
 *
 * G4 -- AIR, Cu, KAPTON, Ti, Al
 * A2 -- SS, NbTi, Epoxy, HeButanol, HeMix
 * */
/* ******************************************** CLASSES *********************************************
 * G4VisAttributes(G4Colour(red, green, blue, opacity))
 * red, green, blue on scale from 0 to 1
 * opacity set to 1 but can be changed
 *
 * G4Tubs("Name", fRMin, fRMax, fDz, fSPhi, fDPhi) -- Tube/Cylinder
 * fRMin = inner radius. If inner radius is 0, tube is filled cylinder.
 * fRMax = outer radius
 * fDz = half length along z-axis, along which it is centered
 * fSPhi = starting phi angle (radians). 0 is +x axis, pi/2 is +y axis.
 * fDPhi = delta angle (radians). 2pi is full cylinder.
 *
 * G4Box("Name", x, y, z) -- 3D Box
 * x, y, z = half-lengths
 *
 * G4LogicalVolume(pSolid, pMaterial, "Name") -- Volume in G4
 * pSolid = G4Tubs or G4Box for example -- shape or volume created
 * pMaterial = material. Found in fNistManager->FindOrBuildMaterial("G4_XXX")
 *
 * G4PVPlacement(pRot, G4ThreeVector(x, y, z), G4LogicalVolume, "Name", pMotherLogical, pMany, pCopyNo) -- Location in Space
 * pRot = volume rotated by pRot relative to pMotherLogical
 * G4ThreeVector = volume trasnlated by it relative to pMotherLogical
 * pMotherLogical = fMyLogic always, I think
 * pMany = false. (It's used for overlapping, but currently not used.)
 * pCopyNo = 1 for everything I can see
 *
 * G4SubtractionSolid("Name", pSolidA, pSolidB)
 * pSolidA, pSolidB = G4Tubs or G4Box, for example -- shape or volume created
 * Subtraction A - B, meaning B is subtracted from A.
 * */
/* ***************************************** OTHER TARGETS ******************************************
 * Cryo target -- G4Sphere("Name", rmin, rmax, phimain, delatphi, thetamin, deltatheta)
 * Options for LH2 or LD2 are imbedded into DetectorSetup.mac and then looked for within Cryo
 *
 * Solid target -- Useful line: if (!fMaterial) give an error message.
 * if fMaterial is something, give it a length (4 options)
 * Cone: G4Cons("Name", innerradius1, outerradius1, innerradius2, outerradius2, halflengthinz, startingphi, deltaphi)
 * Target holder?
 *
 * Structure of the fMaterial concept:
 *
 * G4double trgt_length;
 * if (!fMaterial) {G4cerr<<"A2PolarisedTarget::Construct() Polarised target material not defined. Add in DetectorSetup.mac."<<G4endl;exit(1);}
 * else if (fMaterial == G4NistManager::Instance()->FindOrBuildMaterial("A2_HeButanol")) trgt_length = 20.0*mm;
 * else if (fMaterial == G4NistManager::Instance()->FindOrBuildMaterial("A2_Polarized")) trgt_length = .*mm; // MAYBE, IF IT'S CALLED POLARIZED
 * else {G4cerr<<"A2PolarisedTarget::Construct() Polarised target material not allowed. Change in DetectorSetup.mac."<<G4endl;exit(1);}
 * fLength = trgt_length;
 *
 * In SolidTarget, CELLLogic is made up of the material fMaterial
 * In CryoTarget, LD2ALogic is made up of the material fMaterial
 *
 * */
/* *************************************** ORIGINAL COMMENTS ****************************************
 * // Magnetic field moved from A2DetectorConstruction dglazier
 * static G4bool fieldIsInitialized = false;
 * if(!fieldIsInitialized) {
 * fMagneticField = new A2MagneticField();
 * G4FieldManager* fieldMgr = G4TransportationManager::GetTransportationManager()->GetFieldManager();
 * fieldMgr->SetDetectorField(fMagneticField);
 * fieldMgr->CreateChordFinder(fMagneticField);
 * fieldIsInitialized = true;}
 *
 * There are 2 approximations in the target geometries:
 * When the tubes change thickness, there is no slope joining the different thicknesses.
 * When the ends of the tubes are rounded, this is approximated as a straight 90*degree cylinder.
 * l = length, r = radius, t = thickness.*/

A2PolarizedTarget::A2PolarizedTarget()
{
  fMagneticField = NULL;
}
A2PolarizedTarget::~A2PolarizedTarget()
{
  if(fMagneticField) delete fMagneticField;
}

void A2PolarizedTarget::SetMagneticField(G4String &nameFileFieldMap)
{
  // If nameFileFieldMap is a NULL string then do not set the target magnetic field
  if(nameFileFieldMap.isNull()) {G4cout<<"Warning A2PolarizedTarget::SetMagneticField No field map given, therefore there will be no field!"<<G4endl;return;}
  
  // Create magnetic field
  fMagneticField = new A2MagneticField();
  
  // Read magnetic field map
  // Set this field as default and create trajectory calculator
  // Or, in case of a problem reading the field map, delete fMagneticField
  if(fMagneticField->ReadFieldMap(nameFileFieldMap))
  {
    G4FieldManager* fieldMgr = G4TransportationManager::GetTransportationManager()->GetFieldManager();
    fieldMgr->SetDetectorField(fMagneticField);
    fieldMgr->CreateChordFinder(fMagneticField);
  }
  else
  {
    delete fMagneticField;
  }
}

G4VPhysicalVolume* A2PolarizedTarget::Construct(G4LogicalVolume *MotherLogic) {

 //////////////////////////////////////////////////
 /// Initialization of Mother Volume
 //////////////////////////////////////////////////

 G4bool active = false;
 if (!fMaterial) {G4cerr << "No target material chosen. Please add in DetectorSetup.mac." << G4endl; exit(1);}
 else if (fMaterial == G4NistManager::Instance()->FindOrBuildMaterial("A2_HeButanol")) {G4cout << "A2_HeButanol chosen." << G4endl;}
 else if (fMaterial == G4NistManager::Instance()->FindOrBuildMaterial("G4_POLYSTYRENE")) {G4cout << "Active Target (G4_POLYSTYRENE) chosen." << G4endl; active = true;}
 else {G4cerr << "Target material incorrectly chosen. Please change in DetectorSetup.mac." << G4endl; exit(1);}

 fMotherLogic=MotherLogic;
 fLength=20.0*mm;
 fRadius=0.5*mm;

 G4double l_TRGT = 276.*mm + 67.*mm + 2*mm;
 G4double r_TRGT = 47.97*mm;
  
 // Mother volume, a cylinder with at least a 1-mm space between it and the target in all directions:
 G4Tubs* MyShape=new G4Tubs("TRGT",0.,r_TRGT + 1.0*mm, l_TRGT/2 +1.0*mm,0*deg,360*deg);
 fMyLogic=new G4LogicalVolume(MyShape,fNistManager->FindOrBuildMaterial("G4_AIR"),"TRGT");
 fMyPhysi=new G4PVPlacement(0,G4ThreeVector(0,0, - 20.0*mm/2 - 11.5*mm - 231.5*mm + l_TRGT/2.),fMyLogic,"TRGT",fMotherLogic,false,1);
 fMyLogic->SetVisAttributes (G4VisAttributes::Invisible);

 // Colours with their corresponding materials used in the visualization:
 G4VisAttributes* SSVisAtt= new G4VisAttributes(G4Colour(0.8,0.8,0.8)); // stainless steel (grey)
 G4VisAttributes* CUVisAtt= new G4VisAttributes(G4Colour(0.8,0.6,0.2)); // copper (brown)
 G4VisAttributes* CyanVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,0.0)); // NbTi
 G4VisAttributes* MagentaVisAtt= new G4VisAttributes(G4Colour(1.0,0.0,1.0)); // epoxy
 G4VisAttributes* BlueVisAtt= new G4VisAttributes(G4Colour(0.0,0.0,1.0)); // titanium
 G4VisAttributes* GreenVisAtt= new G4VisAttributes(G4Colour(0.0,1.0,0.0)); // aluminum
 G4VisAttributes* RedVisAtt= new G4VisAttributes(G4Colour(1.0, 0.0, 0.0)); // kapton
 G4VisAttributes* WhiteVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0)); // helium
 G4VisAttributes* PlexiVisAtt= new G4VisAttributes(G4Colour(1.0,1.0, 0.0)); // plexiglass
 G4VisAttributes* BlackVisAtt= new G4VisAttributes(G4Colour(0.0,0.0,0.0)); // black

 //////////////////////////////////////////////////
 /// Active Target
 //////////////////////////////////////////////////

 if (active == true)
 {
     // Plexiglass Tube
     G4double l_PGTube = 140.*mm;
     G4double r_PGTube = 13.*mm; // diameter of 26 mm
     G4double t_PGTube = 2.9*mm; // such that inside is empty
     G4double tubelocation = l_PGTube/2. + 67*mm + 16*mm - l_TRGT/2.; // from old copper cylinder
     G4Tubs* PGTube = new G4Tubs("PGTube", r_PGTube-t_PGTube, r_PGTube, l_PGTube/2., 0*deg, 360*deg);
     G4LogicalVolume* PGTubeLogic = new G4LogicalVolume(PGTube, fNistManager->FindOrBuildMaterial("G4_PLEXIGLASS"), "PGTube");
     new G4PVPlacement(0, G4ThreeVector(0,0,tubelocation),PGTubeLogic,"PGTube",fMyLogic,false,1);
     PGTubeLogic->SetVisAttributes(RedVisAtt);

     // 10 small polystyrene scintillator slices
     G4double l_PSS = 1*mm;
     G4double r_PSS = 1*cm;
     G4double PSS_start = tubelocation+l_PGTube/2.;

     G4Tubs* PSS1 = new G4Tubs("PSS1", 0, r_PSS, l_PSS/2., 0*deg, 360*deg);
     G4LogicalVolume* PSS1Logic = new G4LogicalVolume(PSS1, fNistManager->FindOrBuildMaterial("G4_POLYSTYRENE"), "PSS1");
     new G4PVPlacement(0, G4ThreeVector(0,0,PSS_start),PSS1Logic,"PSS1",fMyLogic,false,1);
     PSS1Logic->SetVisAttributes(MagentaVisAtt);

     G4Tubs* PSS2 = new G4Tubs("PSS2", 0, r_PSS, l_PSS/2., 0*deg, 360*deg);
     G4LogicalVolume* PSS2Logic = new G4LogicalVolume(PSS2, fNistManager->FindOrBuildMaterial("G4_POLYSTYRENE"), "PSS2");
     new G4PVPlacement(0, G4ThreeVector(0,0,PSS_start+2*mm),PSS2Logic,"PSS2",fMyLogic,false,1);
     PSS2Logic->SetVisAttributes(MagentaVisAtt);

     G4Tubs* PSS3 = new G4Tubs("PSS3", 0, r_PSS, l_PSS/2., 0*deg, 360*deg);
     G4LogicalVolume* PSS3Logic = new G4LogicalVolume(PSS3, fNistManager->FindOrBuildMaterial("G4_POLYSTYRENE"), "PSS3");
     new G4PVPlacement(0, G4ThreeVector(0,0,PSS_start+4*mm),PSS3Logic,"PSS3",fMyLogic,false,1);
     PSS3Logic->SetVisAttributes(MagentaVisAtt);

     G4Tubs* PSS4 = new G4Tubs("PSS4", 0, r_PSS, l_PSS/2., 0*deg, 360*deg);
     G4LogicalVolume* PSS4Logic = new G4LogicalVolume(PSS4, fNistManager->FindOrBuildMaterial("G4_POLYSTYRENE"), "PSS4");
     new G4PVPlacement(0, G4ThreeVector(0,0,PSS_start+6*mm),PSS4Logic,"PSS4",fMyLogic,false,1);
     PSS4Logic->SetVisAttributes(MagentaVisAtt);

     G4Tubs* PSS5 = new G4Tubs("PSS5", 0, r_PSS, l_PSS/2., 0*deg, 360*deg);
     G4LogicalVolume* PSS5Logic = new G4LogicalVolume(PSS5, fNistManager->FindOrBuildMaterial("G4_POLYSTYRENE"), "PSS5");
     new G4PVPlacement(0, G4ThreeVector(0,0,PSS_start+8*mm),PSS5Logic,"PSS5",fMyLogic,false,1);
     PSS5Logic->SetVisAttributes(MagentaVisAtt);

     G4Tubs* PSS6 = new G4Tubs("PSS6", 0, r_PSS, l_PSS/2., 0*deg, 360*deg);
     G4LogicalVolume* PSS6Logic = new G4LogicalVolume(PSS6, fNistManager->FindOrBuildMaterial("G4_POLYSTYRENE"), "PSS6");
     new G4PVPlacement(0, G4ThreeVector(0,0,PSS_start+10*mm),PSS6Logic,"PSS6",fMyLogic,false,1);
     PSS6Logic->SetVisAttributes(MagentaVisAtt);

     G4Tubs* PSS7 = new G4Tubs("PSS7", 0, r_PSS, l_PSS/2., 0*deg, 360*deg);
     G4LogicalVolume* PSS7Logic = new G4LogicalVolume(PSS7, fNistManager->FindOrBuildMaterial("G4_POLYSTYRENE"), "PSS7");
     new G4PVPlacement(0, G4ThreeVector(0,0,PSS_start+12*mm),PSS7Logic,"PSS7",fMyLogic,false,1);
     PSS7Logic->SetVisAttributes(MagentaVisAtt);

     G4Tubs* PSS8 = new G4Tubs("PSS8", 0, r_PSS, l_PSS/2., 0*deg, 360*deg);
     G4LogicalVolume* PSS8Logic = new G4LogicalVolume(PSS8, fNistManager->FindOrBuildMaterial("G4_POLYSTYRENE"), "PSS8");
     new G4PVPlacement(0, G4ThreeVector(0,0,PSS_start+14*mm),PSS8Logic,"PSS8",fMyLogic,false,1);
     PSS8Logic->SetVisAttributes(MagentaVisAtt);

     G4Tubs* PSS9 = new G4Tubs("PSS9", 0, r_PSS, l_PSS/2., 0*deg, 360*deg);
     G4LogicalVolume* PSS9Logic = new G4LogicalVolume(PSS9, fNistManager->FindOrBuildMaterial("G4_POLYSTYRENE"), "PSS9");
     new G4PVPlacement(0, G4ThreeVector(0,0,PSS_start+16*mm),PSS9Logic,"PSS9",fMyLogic,false,1);
     PSS9Logic->SetVisAttributes(MagentaVisAtt);

     G4Tubs* PSS10 = new G4Tubs("PSS10", 0, r_PSS, l_PSS/2., 0*deg, 360*deg);
     G4LogicalVolume* PSS10Logic = new G4LogicalVolume(PSS10, fNistManager->FindOrBuildMaterial("G4_POLYSTYRENE"), "PSS10");
     new G4PVPlacement(0, G4ThreeVector(0,0,PSS_start+18*mm),PSS10Logic,"PSS10",fMyLogic,false,1);
     PSS10Logic->SetVisAttributes(MagentaVisAtt);
/*
     // Cone Cap
     G4double l_CONECAP = 2*cm;
     G4double r_CONECAP = 1.25*cm;
     G4Cons* ConeCap = new G4Cons("Cone Cap",0,0.001*cm,0*cm,r_CONECAP,l_CONECAP/2.,0*deg,360*deg);
     G4LogicalVolume* ConeCapLogic = new G4LogicalVolume(ConeCap, fNistManager->FindOrBuildMaterial("G4_PLEXIGLASS"), "ConeCap");
     G4RotationMatrix cone_rotm = G4RotationMatrix();
     cone_rotm.rotateY(180*deg);
     G4ThreeVector conelocation = G4ThreeVector(0,3*cm,tubelocation+l_PGTube/2.+l_CONECAP/2.);
     G4Transform3D cone_transform = G4Transform3D(cone_rotm,conelocation);
     new G4PVPlacement(cone_transform,ConeCapLogic,"ConeCap",fMyLogic,false,1);
     ConeCapLogic->SetVisAttributes(RedVisAtt);

     // Solid Sphere Cap
     G4double r_SPHERECAP = 1.25*cm;
     G4Sphere* SphereCap = new G4Sphere("SphereCap",0,r_SPHERECAP,0*deg,180*deg,0*deg,180*deg);
     G4LogicalVolume* SphereCapLogic = new G4LogicalVolume(SphereCap,fNistManager->FindOrBuildMaterial("G4_PLEXIGLASS"),"SphereCap");
     G4RotationMatrix sphere_rotm = G4RotationMatrix();
     sphere_rotm.rotateX(90*deg);
     G4ThreeVector spherelocation = G4ThreeVector(0,0,tubelocation+l_PGTube/2.);
     G4Transform3D sphere_transform = G4Transform3D(sphere_rotm,spherelocation);
     new G4PVPlacement(sphere_transform,SphereCapLogic,"SphereCap",fMyLogic,false,1);
     SphereCapLogic->SetVisAttributes(RedVisAtt);

     // Hollow Sphere Cap
     G4double t_SPHERECAP = 0.25*cm;
     G4Sphere* HollowSphereCap = new G4Sphere("HollowSphereCap",r_SPHERECAP-t_SPHERECAP,r_SPHERECAP,0*deg,180*deg,0*deg,180*deg);
     G4LogicalVolume* HollowSphereCapLogic = new G4LogicalVolume(HollowSphereCap,fNistManager->FindOrBuildMaterial("G4_PLEXIGLASS"),"HollowSphereCap");
     G4RotationMatrix hollowsphere_rotm = G4RotationMatrix();
     hollowsphere_rotm.rotateX(90*deg);
     G4ThreeVector hollowspherelocation = G4ThreeVector(0,-3*cm,tubelocation+l_PGTube/2.);
     G4Transform3D hollowsphere_transform = G4Transform3D(hollowsphere_rotm,hollowspherelocation);
     new G4PVPlacement(hollowsphere_transform,HollowSphereCapLogic,"HollowSphereCap",fMyLogic,false,1);
     HollowSphereCapLogic->SetVisAttributes(RedVisAtt);

     // Paraboloid Cap
     G4double l_paraboloid = 2*cm;
     G4double paraboloidradius = 1.25*cm;
     G4Paraboloid* ParaboloidCap = new G4Paraboloid("ParaboloidCap",l_paraboloid/2.,0,paraboloidradius);
     G4LogicalVolume* ParaboloidCapLogic = new G4LogicalVolume(ParaboloidCap,fNistManager->FindOrBuildMaterial("G4_PLEXIGLASS"),"ParaboloidCap");
     G4RotationMatrix paraboloid_rotm = G4RotationMatrix();
     paraboloid_rotm.rotateX(180*deg);
     G4ThreeVector paraboloidlocation = G4ThreeVector(0, 6*cm, tubelocation+l_PGTube/2.+l_paraboloid/2.);
     G4Transform3D paraboloid_transform = G4Transform3D(paraboloid_rotm, paraboloidlocation);
     new G4PVPlacement(paraboloid_transform,ParaboloidCapLogic,"ParaboloidCap",fMyLogic,false,1);
     ParaboloidCapLogic->SetVisAttributes(RedVisAtt);

     // Polycone Cap
     G4double numberoflayers = 5;
     const G4double rInner[] = {0*mm,0*mm,0*mm,0*mm,0*mm};
     const G4double rOuter[] = {1.25*cm,1.25*cm,1*cm,0.5*cm,0.1*cm};
     const G4double zslices[] = {0*mm,7*mm,13*mm,19*mm,20*mm};
     G4Polycone* PolyconeCap = new G4Polycone("PolyconeCap",0*deg,360*deg,numberoflayers,zslices,rInner,rOuter);
     G4LogicalVolume* PolyconeCapLogic = new G4LogicalVolume(PolyconeCap,fNistManager->FindOrBuildMaterial("G4_PLEXIGLASS"),"PolyconeCap");
     new G4PVPlacement(0,G4ThreeVector(0,-6*cm,tubelocation+l_PGTube/2.),PolyconeCapLogic,"PolyconeCap",fMyLogic,false,1);
     PolyconeCapLogic->SetVisAttributes(RedVisAtt);
*/
 }

 //////////////////////////////////////////////////
 /// Not Active Polarized Target
 //////////////////////////////////////////////////

 if (active == false)
 {
     //////////////////////////////////////////////////
     /// Inner Cu Cylinders
     //////////////////////////////////////////////////

     // Inner copper cylinder, part A:
     G4double l_CUIA = 10.4*mm;
     G4double r_CUIA = 11.5*mm;
     G4double t_CUIA = 1.5*mm;
     G4Tubs* CUIA=new G4Tubs("CUIA",r_CUIA-t_CUIA,r_CUIA,l_CUIA/2,0*deg,360*deg);
     G4LogicalVolume* CUIALogic=new G4LogicalVolume(CUIA,fNistManager->FindOrBuildMaterial("G4_Cu"),"CUIA");
     new G4PVPlacement(0,G4ThreeVector(0,0,(231.5*mm + 1.0*mm + l_CUIA/2 - l_TRGT/2.)),CUIALogic,"CUIA",fMyLogic,false,1);
     CUIALogic->SetVisAttributes(CUVisAtt);

     // Inner copper cylinder, part B:
     G4double l_CUIB = 1.0*mm;
     G4double r_CUIB = 12.2*mm;
     G4double t_CUIB = 2.2*mm;
     G4Tubs* CUIB=new G4Tubs("CUIB",r_CUIB-t_CUIB,r_CUIB,l_CUIB/2,0*deg,360*deg);
     G4LogicalVolume* CUIBLogic=new G4LogicalVolume(CUIB,fNistManager->FindOrBuildMaterial("G4_Cu"),"CUIB");
     new G4PVPlacement(0,G4ThreeVector(0,0,(231.5*mm + l_CUIB/2 - l_TRGT/2.)),CUIBLogic,"CUIB",fMyLogic,false,1);
     CUIBLogic->SetVisAttributes(CUVisAtt);

     // Inner copper cylinder, part C:
     G4double l_CUIC = 9.0*mm;
     G4double r_CUIC = 12.2*mm;
     G4double t_CUIC = 0.5*mm;
     G4Tubs* CUIC=new G4Tubs("CUIC",r_CUIC-t_CUIC,r_CUIC,l_CUIC/2,0*deg,360*deg);
     G4LogicalVolume* CUICLogic=new G4LogicalVolume(CUIC,fNistManager->FindOrBuildMaterial("G4_Cu"),"CUIC");
     new G4PVPlacement(0,G4ThreeVector(0,0,(231.5*mm - l_CUIC/2 - l_TRGT/2.)),CUICLogic,"CUIC",fMyLogic,false,1);
     CUICLogic->SetVisAttributes(CUVisAtt);

     // Inner copper cylinder, part D:
     G4double l_CUID = 215.*mm - 67*mm;
     G4double r_CUID = 12.5*mm;
     G4double t_CUID = 0.3*mm;
     G4Tubs* CUID=new G4Tubs("CUID",r_CUID-t_CUID,r_CUID,l_CUID/2,0*deg,360*deg);
     G4LogicalVolume* CUIDLogic=new G4LogicalVolume(CUID,fNistManager->FindOrBuildMaterial("G4_Cu"),"CUID");
     new G4PVPlacement(0,G4ThreeVector(0,0,(l_CUID/2 + 67*mm + 16*mm - l_TRGT/2.)),CUIDLogic,"CUID",fMyLogic,false,1);
     CUIDLogic->SetVisAttributes(CUVisAtt);

     //////////////////////////////////////////////////
     /// Butanol Target and Kapton Cell
     //////////////////////////////////////////////////

     // The butanol target is encased in kapton. This cup has holes in it, which are ingnored.
     // Kapton outside cell, part A:
     G4double l_KAPA = 0.6*mm;
     G4double r_KAPA = 10.505*mm;
     G4Tubs* KAPA=new G4Tubs("KAPA",0,r_KAPA,l_KAPA/2,0*deg,360*deg);
     G4LogicalVolume* KAPALogic=new G4LogicalVolume(KAPA,fNistManager->FindOrBuildMaterial("G4_KAPTON"),"KAPA");
     new G4PVPlacement(0,G4ThreeVector(0,0,(l_KAPA/2 + 20.0*mm + 11.5*mm + 231.5*mm - l_TRGT/2.)),KAPALogic,"KAPA",fMyLogic,false,1);
     KAPALogic->SetVisAttributes(RedVisAtt);

     // Kapton outside cell, part B:
     G4double l_KAPB = 19.0*mm;
     G4double r_KAPB = 10.505*mm;
     G4double t_KAPB = 0.6*mm;
     G4Tubs* KAPB=new G4Tubs("KAPB",r_KAPB-t_KAPB,r_KAPB,l_KAPB/2,0*deg,360*deg);
     G4LogicalVolume* KAPBLogic=new G4LogicalVolume(KAPB,fNistManager->FindOrBuildMaterial("G4_KAPTON"),"KAPB");
     new G4PVPlacement(0,G4ThreeVector(0,0,(l_KAPB/2 + 1.0*mm + 11.5*mm + 231.5*mm - l_TRGT/2.)),KAPBLogic,"KAPB",fMyLogic,false,1);
     KAPBLogic->SetVisAttributes(RedVisAtt);

     // Kapton outside cell, part C:
     G4double l_KAPC = 1.0*mm;
     G4double r_KAPC = 12.43*mm;
     G4double t_KAPC = 2.5*mm;
     G4Tubs* KAPC=new G4Tubs("KAPC",r_KAPC-t_KAPC,r_KAPC,l_KAPC/2,0*deg,360*deg);
     G4LogicalVolume* KAPCLogic=new G4LogicalVolume(KAPC,fNistManager->FindOrBuildMaterial("G4_KAPTON"),"KAPC");
     new G4PVPlacement(0,G4ThreeVector(0,0,(l_KAPC/2 + 11.5*mm + 231.5*mm - l_TRGT/2.)),KAPCLogic,"KAPC",fMyLogic,false,1);
     KAPCLogic->SetVisAttributes(RedVisAtt);

     // Kapton outside cell, part D:
     G4double l_KAPD = 10.0*mm;
     G4double r_KAPD = 12.43*mm;
     G4double t_KAPD = 0.93*mm;
     G4Tubs* KAPD=new G4Tubs("KAPD",r_KAPD-t_KAPD,r_KAPD,l_KAPD/2,0*deg,360*deg);
     G4LogicalVolume* KAPDLogic=new G4LogicalVolume(KAPD,fNistManager->FindOrBuildMaterial("G4_KAPTON"),"KAPD");
     new G4PVPlacement(0,G4ThreeVector(0,0,(l_KAPD/2 + 1.5*mm + 231.5*mm - l_TRGT/2.)),KAPDLogic,"KAPD",fMyLogic,false,1);
     KAPDLogic->SetVisAttributes(RedVisAtt);

     // Target cell, butanol, 60% filling:
     G4double l_BTRGT = 20.0*mm;
     G4double r_BTRGT = 9.905*mm;
     G4Tubs* BTRGT=new G4Tubs("BTRGT",0,r_BTRGT,l_BTRGT/2,0*deg,360*deg);
     G4LogicalVolume* BTRGTLogic=new G4LogicalVolume(BTRGT,fNistManager->FindOrBuildMaterial("A2_HeButanol"),"BTRGT"); // changed from A2_HeButanol to fMaterial
     //G4LogicalVolume* BTRGTLogic=new G4LogicalVolume(BTRGT,fNistManager->FindOrBuildMaterial("A2_HeMix"),"BTRGT");
     new G4PVPlacement(0,G4ThreeVector(0,0,(l_BTRGT/2 + 11.5*mm + 231.5*mm - l_TRGT/2.)),BTRGTLogic,"BTRGT",fMyLogic,false,1);
     BTRGTLogic->SetVisAttributes(MagentaVisAtt);

     //////////////////////////////////////////////////
     /// Helium between Inner SS and Cu Cylinders
     //////////////////////////////////////////////////

     // He, part A:
     G4double l_HEA = 3.1*mm;
     G4double r_HEA = 19.5*mm;
     G4Tubs* HEA=new G4Tubs("HEA",0,r_HEA,l_HEA/2,0*deg,360*deg);
     G4LogicalVolume* HEALogic=new G4LogicalVolume(HEA,fNistManager->FindOrBuildMaterial("A2_HeMix"),"HEA");
     new G4PVPlacement(0,G4ThreeVector(0,0,(67*mm + 200*mm - 0.3*mm - l_HEA/2 - l_TRGT/2.)),HEALogic,"HEA",fMyLogic,false,1);
     HEALogic->SetVisAttributes(WhiteVisAtt);

     // He, part B:
     G4double l_HEB = 19.575*mm;
     G4double r_HEB = 19.5*mm;
     G4double t_HEB = 8.995*mm;
     G4Tubs* HEB=new G4Tubs("HEB",r_HEB-t_HEB,r_HEB,l_HEB/2,0*deg,360*deg);
     G4LogicalVolume* HEBLogic=new G4LogicalVolume(HEB,fNistManager->FindOrBuildMaterial("A2_HeMix"),"HEB");
     new G4PVPlacement(0,G4ThreeVector(0,0,(l_HEB/2 + 1*mm + 11.5*mm + 231.5*mm - l_TRGT/2.)),HEBLogic,"HEB",fMyLogic,false,1);
     HEBLogic->SetVisAttributes(WhiteVisAtt);

     // He, part C:
     G4double l_HEC = 11.0*mm;
     G4double r_HEC = 19.5*mm;
     G4double t_HEC = 7.07*mm;
     G4Tubs* HEC=new G4Tubs("HEC",r_HEC-t_HEC,r_HEC,l_HEC/2,0*deg,360*deg);
     G4LogicalVolume* HECLogic=new G4LogicalVolume(HEC,fNistManager->FindOrBuildMaterial("A2_HeMix"),"HEC");
     new G4PVPlacement(0,G4ThreeVector(0,0,(231.5*mm + 1.5*mm + l_HEC/2 - l_TRGT/2.)),HECLogic,"HEC",fMyLogic,false,1);
     HECLogic->SetVisAttributes(WhiteVisAtt);

     // He, part D:
     G4double l_HED = 0.5*mm;
     G4double r_HED = 19.5*mm;
     G4double t_HED = 8.0*mm;
     G4Tubs* HED=new G4Tubs("HED",r_HED-t_HED,r_HED,l_HED/2,0*deg,360*deg);
     G4LogicalVolume* HEDLogic=new G4LogicalVolume(HED,fNistManager->FindOrBuildMaterial("A2_HeMix"),"HED");
     new G4PVPlacement(0,G4ThreeVector(0,0,(231.5*mm + 1*mm + l_HED/2 - l_TRGT/2.)),HEDLogic,"HED",fMyLogic,false,1);
     HEDLogic->SetVisAttributes(WhiteVisAtt);

     // He, part E:
     G4double l_HEE = 1.0*mm;
     G4double r_HEE = 19.5*mm;
     G4double t_HEE = 7.3*mm;
     G4Tubs* HEE=new G4Tubs("HEE",r_HEE-t_HEE,r_HEE,l_HEE/2,0*deg,360*deg);
     G4LogicalVolume* HEELogic=new G4LogicalVolume(HEE,fNistManager->FindOrBuildMaterial("A2_HeMix"),"HEE");
     new G4PVPlacement(0,G4ThreeVector(0,0,(231.5*mm + l_HEE/2 - l_TRGT/2.)),HEELogic,"HEE",fMyLogic,false,1);
     HEELogic->SetVisAttributes(WhiteVisAtt);

     // He, part F:
     G4double l_HEF = 122.0*mm;
     G4double r_HEF = 19.5*mm;
     G4double t_HEF = 7.0*mm;
     G4Tubs* HEF=new G4Tubs("HEF",r_HEF-t_HEF,r_HEF,l_HEF/2,0*deg,360*deg);
     G4LogicalVolume* HEFLogic=new G4LogicalVolume(HEF,fNistManager->FindOrBuildMaterial("A2_HeMix"),"HEF");
     new G4PVPlacement(0,G4ThreeVector(0,0,(231.5*mm - l_HEF/2 - l_TRGT/2.)),HEFLogic,"HEF",fMyLogic,false,1);
     HEFLogic->SetVisAttributes(WhiteVisAtt);
 }

 //////////////////////////////////////////////////
 /// Cylinders, from the outside working in
 //////////////////////////////////////////////////

 // Outer stainless steel cylinder:
 G4double l_SSO = 233.5*mm;
 G4double r_SSO = 32.5*mm;
 G4double t_SSO = 0.5*mm;
 G4Tubs* SSO=new G4Tubs("SSO",r_SSO-t_SSO,r_SSO,l_SSO/2,0*deg,360*deg);
 G4LogicalVolume* SSOLogic=new G4LogicalVolume(SSO,fNistManager->FindOrBuildMaterial("A2_SS"),"SSO");
 new G4PVPlacement(0,G4ThreeVector(0,0,(l_SSO/2 + 42.5*mm + 67*mm - l_TRGT/2.)),SSOLogic,"SSO",fMyLogic,false,1);
 SSOLogic->SetVisAttributes(SSVisAtt);

 // Outer copper cylinder:
 G4double l_CUO = 222.5*mm;
 G4double r_CUO = 29.0*mm;
 G4double t_CUO = 0.5*mm;
 G4Tubs* CUO=new G4Tubs("CUO",r_CUO-t_CUO,r_CUO,l_CUO/2,0*deg,360*deg);
 G4LogicalVolume* CUOLogic=new G4LogicalVolume(CUO,fNistManager->FindOrBuildMaterial("G4_Cu"),"CUO");
 new G4PVPlacement(0,G4ThreeVector(0,0,(l_CUO/2 + 42.5*mm + 67*mm - l_TRGT/2.)),CUOLogic,"CUO",fMyLogic,false,1);
 CUOLogic->SetVisAttributes(CUVisAtt);

 // The coils are wrapped around part of a middle copper cylinder, which has 3 different thicknesses, 0.3, 0.7, and 1.2, labelled A, B, and C respectively, going up-beam.
 // Middle copper cylinder, part A, which the coils are wrapped around:
 G4double l_CUMA = 136.0*mm;
 G4double r_CUMA = 23.6*mm;
 G4double t_CUMA = 0.3*mm;
 G4Tubs* CUMA=new G4Tubs("CUMA",r_CUMA-t_CUMA,r_CUMA,l_CUMA/2,0*deg,360*deg);
 G4LogicalVolume* CUMALogic=new G4LogicalVolume(CUMA,fNistManager->FindOrBuildMaterial("G4_Cu"),"CUMA");
 new G4PVPlacement(0,G4ThreeVector(0,0,(l_CUMA/2 + 119*mm + 67*mm - l_TRGT/2.)),CUMALogic,"CUMA",fMyLogic,false,1);
 CUMALogic->SetVisAttributes(CUVisAtt);

 // Middle copper cylinder, part B:
 G4double l_CUMB = 33.5*mm;
 G4double r_CUMB = 24.0*mm;
 G4double t_CUMB = 0.7*mm;
 G4Tubs* CUMB=new G4Tubs("CUMB",r_CUMB-t_CUMB,r_CUMB,l_CUMB/2,0*deg,360*deg);
 G4LogicalVolume* CUMBLogic=new G4LogicalVolume(CUMB,fNistManager->FindOrBuildMaterial("G4_Cu"),"CUMB");
 new G4PVPlacement(0,G4ThreeVector(0,0,(l_CUMB/2 + 85.5*mm + 67*mm - l_TRGT/2.)),CUMBLogic,"CUMB",fMyLogic,false,1);
 CUMBLogic->SetVisAttributes(CUVisAtt);

 // Middle copper cylinder, part C:
 G4double l_CUMC = 57.7*mm;
 G4double r_CUMC = 24.5*mm;
 G4double t_CUMC = 1.2*mm;
 G4Tubs* CUMC=new G4Tubs("CUMC",r_CUMC-t_CUMC,r_CUMC,l_CUMC/2,0*deg,360*deg);
 G4LogicalVolume* CUMCLogic=new G4LogicalVolume(CUMC,fNistManager->FindOrBuildMaterial("G4_Cu"),"CUMC");
 new G4PVPlacement(0,G4ThreeVector(0,0,(l_CUMC/2 + 27.8*mm + 67*mm - l_TRGT/2.)),CUMCLogic,"CUMC",fMyLogic,false,1);
 CUMCLogic->SetVisAttributes(CUVisAtt);

 // Inner stainless steel cylinder, going up-beam, part A:
 G4double l_SSIA = 122.0*mm;
 G4double r_SSIA = 19.8*mm;
 G4double t_SSIA = 0.3*mm;
 G4Tubs* SSIA=new G4Tubs("SSIA",r_SSIA-t_SSIA,r_SSIA,l_SSIA/2,0*deg,360*deg);
 G4LogicalVolume* SSIALogic=new G4LogicalVolume(SSIA,fNistManager->FindOrBuildMaterial("A2_SS"),"SSIA");
 new G4PVPlacement(0,G4ThreeVector(0,0,(l_SSIA/2 + 78*mm + 67*mm - l_TRGT/2.)),SSIALogic,"SSIA",fMyLogic,false,1);
 SSIALogic->SetVisAttributes(SSVisAtt);

 // Inner stainless steel cylinder, part B:
 G4double l_SSIB = 56.75*mm;
 G4double r_SSIB = 20.0*mm;
 G4double t_SSIB = 0.5*mm;
 G4Tubs* SSIB=new G4Tubs("SSIB",r_SSIB-t_SSIB,r_SSIB,l_SSIB/2,0*deg,360*deg);
 G4LogicalVolume* SSIBLogic=new G4LogicalVolume(SSIB,fNistManager->FindOrBuildMaterial("A2_SS"),"SSIB");
 new G4PVPlacement(0,G4ThreeVector(0,0,(l_SSIB/2 + 21.25*mm + 67*mm - l_TRGT/2.)),SSIBLogic,"SSIB",fMyLogic,false,1);
 SSIBLogic->SetVisAttributes(SSVisAtt);

 //////////////////////////////////////////////////
 /// Windows, from outside going up-beam
 //////////////////////////////////////////////////

 // The cylinders with the Ti windows attached were rounded on the down-beam end. This makes the approximation of 90*deg corners.
 // Outer Ti window:
 G4double r_TIOW = 12.5*mm;
 G4double t_TIOW = 0.02*mm;
 G4Tubs* TIOW=new G4Tubs("TIOW",0,r_TIOW,t_TIOW/2,0*deg,360*deg);
 G4LogicalVolume* TIOWLogic=new G4LogicalVolume(TIOW,fNistManager->FindOrBuildMaterial("G4_Ti"),"TIOW");
 new G4PVPlacement(0,G4ThreeVector(0,0,(t_TIOW/2 + 276*mm + 67*mm - l_TRGT/2.)),TIOWLogic,"TIOW",fMyLogic,false,1);
 TIOWLogic->SetVisAttributes(BlueVisAtt);

 // Copper ring around outer Ti window:
 G4double r_CUOW = 32.5*mm;
 G4double t_CUOW = 0.40*mm;
 G4Tubs* CUOW=new G4Tubs("CUOW",r_TIOW,r_CUOW,t_CUOW/2,0*deg,360*deg);
 G4LogicalVolume* CUOWLogic=new G4LogicalVolume(CUOW,fNistManager->FindOrBuildMaterial("G4_Cu"),"CUOW");
 new G4PVPlacement(0,G4ThreeVector(0,0,(t_CUOW/2 + 276*mm + 67*mm - l_TRGT/2.)),CUOWLogic,"CUOW",fMyLogic,false,1);
 CUOWLogic->SetVisAttributes(CUVisAtt);

 // Outer Al window:
 G4double r_ALOW = 29.0*mm;
 G4double t_ALOW = 0.01*mm;
 G4Tubs* ALOW=new G4Tubs("ALOW",0,r_ALOW,t_ALOW/2,0*deg,360*deg);
 G4LogicalVolume* ALOWLogic=new G4LogicalVolume(ALOW,fNistManager->FindOrBuildMaterial("G4_Al"),"ALOW");
 new G4PVPlacement(0,G4ThreeVector(0,0,(t_ALOW/2 + 265*mm + 67*mm - l_TRGT/2.)),ALOWLogic,"ALOW",fMyLogic,false,1);
 ALOWLogic->SetVisAttributes(GreenVisAtt);

 // Inner Al window:
 G4double r_ALIW = 23.7*mm;
 G4double t_ALIW = 0.01*mm;
 G4Tubs* ALIW=new G4Tubs("ALIW",0,r_ALIW,t_ALIW/2,0*deg,360*deg);
 G4LogicalVolume* ALIWLogic=new G4LogicalVolume(ALIW,fNistManager->FindOrBuildMaterial("G4_Al"),"ALIW");
 new G4PVPlacement(0,G4ThreeVector(0,0,(t_ALIW/2 + 256.5*mm + 67*mm - l_TRGT/2.)),ALIWLogic,"ALIW",fMyLogic,false,1);
 ALIWLogic->SetVisAttributes(GreenVisAtt);

 // Small copper bit between coils and inner Al window, Part A:
 G4double l_CUBA = 0.5*mm;
 G4double r_CUBA = 23.7*mm;
 G4double t_CUBA = 1.2*mm;
 G4Tubs* CUBA=new G4Tubs("CUBA",r_CUBA-t_CUBA,r_CUBA,l_CUBA/2,0*deg,360*deg);
 G4LogicalVolume* CUBALogic=new G4LogicalVolume(CUBA,fNistManager->FindOrBuildMaterial("G4_Cu"),"CUBA");
 new G4PVPlacement(0,G4ThreeVector(0,0,(l_CUBA/2 + 256.5*mm + 67*mm - l_TRGT/2.)),CUBALogic,"CUBA",fMyLogic,false,1);
 CUBALogic->SetVisAttributes(CUVisAtt);

 // Small copper bit between coils and inner Al window, Part B:
 G4double l_CUBB = 3.18*mm;
 G4double r_CUBB = 23.0*mm;
 G4double t_CUBB = 0.5*mm;
 G4Tubs* CUBB=new G4Tubs("CUBB",r_CUBB-t_CUBB,r_CUBB,l_CUBB/2,0*deg,360*deg);
 G4LogicalVolume* CUBBLogic=new G4LogicalVolume(CUBB,fNistManager->FindOrBuildMaterial("G4_Cu"),"CUBB");
 new G4PVPlacement(0,G4ThreeVector(0,0,(256.5*mm + 67*mm - l_CUBB/2 - l_TRGT/2.)),CUBBLogic,"CUBB",fMyLogic,false,1);
 CUBBLogic->SetVisAttributes(CUVisAtt);

 // Middle Ti window:
 G4double r_TIMW = 11.0*mm;
 G4double t_TIMW = 0.02*mm;
 G4Tubs* TIMW=new G4Tubs("TIMW",0,r_TIMW,t_TIMW/2,0*deg,360*deg);
 G4LogicalVolume* TIMWLogic=new G4LogicalVolume(TIMW,fNistManager->FindOrBuildMaterial("G4_Ti"),"TIMW");
 new G4PVPlacement(0,G4ThreeVector(0,0,(-t_TIMW/2 + 267*mm - l_TRGT/2.)),TIMWLogic,"TIMW",fMyLogic,false,1);
 TIMWLogic->SetVisAttributes(BlueVisAtt);

 // Stainless steel ring around middle Ti window:
 G4double r_SSIW = 19.8*mm;
 G4double t_SSIW = 0.30*mm;
 G4Tubs* SSIW=new G4Tubs("SSIW",r_TIMW,r_SSIW,t_SSIW/2,0*deg,360*deg);
 G4LogicalVolume* SSIWLogic=new G4LogicalVolume(SSIW,fNistManager->FindOrBuildMaterial("A2_SS"),"SSIW");
 new G4PVPlacement(0,G4ThreeVector(0,0,(-t_SSIW/2 + 267*mm - l_TRGT/2.)),SSIWLogic,"SSIW",fMyLogic,false,1);
 SSIWLogic->SetVisAttributes(SSVisAtt);

 // Inner Ti window:
 G4double r_TIIW = 10.0*mm;
 G4double t_TIIW = 0.02*mm;
 G4Tubs* TIIW=new G4Tubs("TIIW",0,r_TIIW,t_TIIW/2,0*deg,360*deg);
 G4LogicalVolume* TIIWLogic=new G4LogicalVolume(TIIW,fNistManager->FindOrBuildMaterial("G4_Ti"),"TIIW");
 new G4PVPlacement(0,G4ThreeVector(0,0,(t_TIIW/2 + 11.5*mm + 231.5*mm - l_TRGT/2.)),TIIWLogic,"TIIW",fMyLogic,false,1);
 TIIWLogic->SetVisAttributes(BlueVisAtt);

 //////////////////////////////////////////////////
 /// Magnetic Field: Solenoidal or Saddle
 //////////////////////////////////////////////////

 if (fTypeMagneticCoils == G4String("Solenoidal") || fTypeMagneticCoils == G4String("solenoidal")) {
     // The coils are approximated to be 3 individual layers of material, 460 microm NiTi, 340 microm Cu, 200 microm epoxy.
     // NbTi coil layer:
     G4double l_NbTiC = 136.0*mm;
     G4double r_NbTiC = 24.6*mm;
     G4double t_NbTiC = 0.460*mm;
     G4Tubs* NbTiC=new G4Tubs("NbTiC",r_NbTiC-t_NbTiC,r_NbTiC,l_NbTiC/2,0*deg,360*deg);
     G4LogicalVolume* NbTiCLogic=new G4LogicalVolume(NbTiC,fNistManager->FindOrBuildMaterial("A2_NbTi"),"NbTiC");
     new G4PVPlacement(0,G4ThreeVector(0,0,(l_NbTiC/2 + 119*mm + 67*mm - l_TRGT/2.)),NbTiCLogic,"NbTiC",fMyLogic,false,1);
     NbTiCLogic->SetVisAttributes(CyanVisAtt);

     // Copper coil layer:
     G4double l_CUC = 136.0*mm;
     G4double r_CUC = 24.14*mm;
     G4double t_CUC = 0.340*mm;
     G4Tubs* CUC=new G4Tubs("CUC",r_CUC-t_CUC,r_CUC,l_CUC/2,0*deg,360*deg);
     G4LogicalVolume* CUCLogic=new G4LogicalVolume(CUC,fNistManager->FindOrBuildMaterial("G4_Cu"),"CUC");
     new G4PVPlacement(0,G4ThreeVector(0,0,(l_CUC/2 + 119*mm + 67*mm - l_TRGT/2.)),CUCLogic,"CUC",fMyLogic,false,1);
     CUCLogic->SetVisAttributes(CUVisAtt);

     // Epoxy coil layer:
     G4double l_EPC = 136.0*mm;
     G4double r_EPC = 23.8*mm;
     G4double t_EPC = 0.200*mm;
     G4Tubs* EPC=new G4Tubs("EPC",r_EPC-t_EPC,r_EPC,l_EPC/2,0*deg,360*deg);
     G4LogicalVolume* EPCLogic=new G4LogicalVolume(EPC,fNistManager->FindOrBuildMaterial("A2_Epoxy"),"EPC");
     new G4PVPlacement(0,G4ThreeVector(0,0,(l_EPC/2 + 119*mm + 67*mm - l_TRGT/2.)),EPCLogic,"EPC",fMyLogic,false,1);
     EPCLogic->SetVisAttributes(MagentaVisAtt);
 }

 else if (fTypeMagneticCoils == G4String("Saddle") || fTypeMagneticCoils == G4String("saddle")) {
     // Box-splitter
     G4double xBoxSplitter = 25*mm;
     G4double yBoxSplitter = 0.1*mm;
     G4double zBoxSplitter = 136*mm;
     G4Box *boxSplitter = new G4Box("boxSplitter", xBoxSplitter, yBoxSplitter, zBoxSplitter);

     // Layer1
     G4double lTube1 = 135.48*mm;
     G4double rTube1 =  24.46*mm;
     G4double tTube1 =   0.46*mm;
     G4Tubs *tube1 = new G4Tubs("tube1", rTube1 - tTube1, rTube1, lTube1/2, 0*deg, 360*deg);
     G4double xBox1 = 6*mm;
     G4double yBox1 = 50*mm;
     G4double zBox1 = 36*mm;
     G4Box *box1 = new G4Box("box1", xBox1, yBox1, zBox1);
     G4SubtractionSolid *tube_box1 = new G4SubtractionSolid("tube1-box1", tube1, box1);
     G4SubtractionSolid *layer1 = new G4SubtractionSolid("layer1", tube_box1, boxSplitter);
     G4LogicalVolume* logicSaddleCoilsLayer1 = new G4LogicalVolume(layer1, fNistManager->FindOrBuildMaterial("A2_NbTi"), "logicSaddleCoilsLayer1");
     new G4PVPlacement(0, G4ThreeVector(0.,0.,(lTube1/2. + 119*mm + 67*mm - l_TRGT/2.)), logicSaddleCoilsLayer1, "physSaddleCoilsLayer1", fMyLogic, false, 1);
     logicSaddleCoilsLayer1->SetVisAttributes(CyanVisAtt);

     // Layer2
     G4double lTube2 = 135.88*mm;
     G4double rTube2 =  24.92*mm;
     G4double tTube2 =   0.46*mm;
     G4Tubs *tube2 = new G4Tubs("tube2", rTube2 - tTube2, rTube2, lTube2/2, 0*deg, 360*deg);
     G4double xBox2 = 20*mm;
     G4double yBox2 = 50*mm;
     G4double zBox2 = 50*mm;
     G4Box *box2 = new G4Box("box2", xBox2, yBox2, zBox2);
     G4SubtractionSolid *tube_box2 = new G4SubtractionSolid("tube2-box2", tube2, box2);
     G4SubtractionSolid *layer2 = new G4SubtractionSolid("layer2", tube_box2, boxSplitter);
     G4LogicalVolume* logicSaddleCoilsLayer2 = new G4LogicalVolume(layer2, fNistManager->FindOrBuildMaterial("A2_NbTi"), "logicSaddleCoilsLayer2");
     new G4PVPlacement(0, G4ThreeVector(0.,0.,(lTube2/2 + 119*mm + 67*mm - l_TRGT/2.)), logicSaddleCoilsLayer2, "physSaddleCoilsLayer2", fMyLogic, false, 1);
     logicSaddleCoilsLayer2->SetVisAttributes(CyanVisAtt);
 }
 return fMyPhysi;
}
