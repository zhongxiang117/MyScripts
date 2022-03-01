#!/usr/bin/perl -w

# version 0.1.0     : Mar 1, 2022

use strict;
use MdmDiscoveryScript;
use Time::HiRes;

# Check that script is being run from within Discovery Studio
if ( DiscoveryScript::IsStandalone() )
{
    print "Error: This example script is intended to run inside Discovery Studio.\n";
    print "Error: Please load Discovery Studio and try again.\n";
    exit;
}

# Get the last active MDM document (Molecule Window) in Discovery Studio.
my $document = DiscoveryScript::LastActiveDocument(MdmModelType);
if (!$document)
{
    printf "Error: cannot find the correct window\n";
    exit;
}

my $allmols = $document->AllMolecules;
my $nummols = $allmols->Count;


printf "Note: Number of Molecules: %s\n", $nummols;
$document->HideAll();
$document->Center();

for ( my $i = 0 ; $i < $nummols; $i++ )
{
    my $mol = $allmols->Item($i);
    $mol->SetVisible();
    #sleep(1); #sleep on seconds only, int
    Time::HiRes::sleep(1); #unit seconds
    #Time::HiRes::usleep(1);  #unit microsecond.
    $mol->Offline();
}


#print "Note: Done: Closing the document.\n";
#$document->Close();