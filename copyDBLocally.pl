use lib '/home/grigoryanlab/library/ProCEDe/branches/gevorg-branch/modules/';
use GENERAL;
use Getopt::Long;
use File::Path;

my %opts = getInput();
umask 0002;                 # make sure all output is readible/writeable by group and owner
my $grp = "grigoryanlab";   # group name under which everything will be created

# create a local path if needed
my $loc = "/data/scratch/grigoryanlab/localDBs";
if (! -d $loc) {
  File::Path::make_path($loc, {group => $grp});
  GENERAL::csystem("chmod g+s $loc");                    # set the sticky bit for the group
}

# set up flags
my $sync = 0;             # does the database need to be synced?
$sync = 1 if ($opts{f});
my $pathDB = "$loc/$opts{n}";

# get file lock (enter exclusive part of the code)
my $lockFile = "$loc/.$opts{n}.lock";
my $lock = GENERAL::getFileLock($lockFile);

# check if progress of a previous copy was aborted -- if yes, need to sync
my $progFile = "$loc/.$opts{n}.progress";
if (-e $progFile) {
  $sync = 1;
  $opts{c} = 1;  # clear the old copy just in case, since something wrong happened last time
} else {
  GENERAL::csystem("touch $progFile");
}

# does the DB already exist? If no, need to sync
if (!$sync && (! -e $pathDB)) {
  $sync = 1;
}

# sync if need to
if ($sync) {
  # clear old version if asked
  if ((-e $pathDB) && ($opts{c})) {
    GENERAL::csystem("rm -r $pathDB");
  }

  # create directory if it is not there already
  File::Path::make_path($pathDB, {group => $grp}) if (! -e $pathDB);

  # split the network list file into the base name and file names, and create the new list
  my ($nameList, $base) = splitList($opts{l}, $pathDB);

  # local (or NFS) copy or through anthill?
  my $pref = (-d $base) ? "" : "anthill:";

  # sync
  GENERAL::csystem("rsync -rlt --no-p --no-g --chmod=Dg+x,g+w,g+r,Fg+r,g+w --files-from=$nameList $pref$base $pathDB");
  GENERAL::crm($nameList);
}

# print resulting list
printf("$pathDB/list\n");

# finish up -- remove progress file and remove lock
GENERAL::crm($progFile);
GENERAL::releaseFileLock($lock);

sub getInput {
  my $usage = GENERAL::usage("Copies a MaDCaT/MASTER database locally. Safe to call simultaneously from many jobs.",
            "-n", "database name",
            "Optional", "",
            "-l", "if the database name is new (unrecognized), specify a list file with paths to the network location",
            "-f", "force copying even if the database already exists (for updates, for example)",
            "-c", "if used with -f, will clear out the database first if it exists");

  my %opts;
  GetOptions (\%opts, "n=s", "l=s", "f", "c");

  if (!defined($opts{n})) {
    die $usage;
  }

  my %knownDBs = ("test", "/u/gevorg/work/tmp/test.list", "bc-30", "/home/grigoryanlab/library/MASTER/lists/bc-30.list",
                  "bc-30-sc-20140815", "/home/ironfs/scratch/grigoryanlab/fzheng/designScore/bc-30-sc-20140815/list");
  if (defined($knownDBs{$opts{n}})) {
    $opts{l} = $knownDBs{$opts{n}};
  } else {
    GENERAL::assert(defined($opts{l}), "'$opts{n}' is not a recognized database, so -l must be specified");
  }

  return %opts;
}

sub getLocalPath {
  return $loc;
}

sub splitList {
  my $olist = shift;
  my $pathDB = shift;
  my $nlist = "$pathDB/list";
  my $nameList = $nlist.".names";
  my $base = undef;

  my $ifh = GENERAL::GetInFH($olist);
  my $ofh = GENERAL::GetOutFH($nlist);
  my $ofh1 = GENERAL::GetOutFH($nameList);
  foreach my $line (<$ifh>) {
    $line = GENERAL::Trim($line);
    next if ($line eq "");
    GENERAL::assert(($line =~ /^(.+)\/([^\/]+\/[^\/]+)$/) ? 1 : 0, "could not parse base/file names out of line '$line'");
    $base = $1 if (!defined($base));
    GENERAL::assert($base eq $1, "difference base paths in files of list '$olist'");
    $ofh->printf("$pathDB/$2\n");
    $ofh1->printf("$2\n");
  }
  close($ofh1);
  close($ofh);
  close($ifh);

  return ($nameList, $base);
}

sub changePermissions {
  my $file = shift;
}

sub getGroup {
  my $file = shift;
  return getgrgid((stat($file))[5]);
}
