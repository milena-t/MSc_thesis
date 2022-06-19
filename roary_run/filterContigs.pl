#!/usr/bin/perl
#This program was created by Bo Segerman, National veterinary Institute, Uppsala, Sweden, bo.segerman@sva.se


############## ARGUMENT HANDLING #################
$treshold= $ARGV[0]; # treshold
$arg_path= $ARGV[1];  #file

$overwrite= lc($ARGV[2]);  #overwrite optional pass 'overwrite' as third argument


$filecount = "";

###################################################

if( -f $arg_path){
 &runFiltering($arg_path); 
}
elsif( -d $arg_path){
  opendir(DIR, $arg_path) or die "Cannot open $arg_path \n";
  
  if(substr($arg_path,-1) ne "/"){
    $arg_path = $arg_path . "/";
  }
  $filecount =0;
  while($file = readdir(DIR)){ 
  


    if( -f $arg_path . $file && ( substr($file,-3) eq ".fa" || substr($file,-6) eq ".fa.gz" || substr($file,-4) eq ".fna" || substr($file,-7) eq ".fna.gz" || substr($file,-6) eq ".fasta" || substr($file,-9) eq ".fasta.gz" )){
       $filecount++;
       &runFiltering($arg_path . $file); 
    }
  }
  closedir( DIR);

}
else{
  die("Not a valid file or directory.. Usage: filterContigs.pl cutofsize file_or_dir_path  optional('overwrite')");
}

###################################################


sub runFiltering(){

$path = shift(@_);

############## FILES #################

if( lc(substr($path,-3)) eq ".gz"){
  open(HANDLE_INFILE,"zcat $path |") || die("FATAL ERROR: coluld not open file $path");
}
else{
  open(HANDLE_INFILE,$path) || die("FATAL ERROR: coluld not open file $path");
}

$outpath = &addtagtopath($path,".filter");

if( lc(substr($outpath,-3)) eq ".gz"){
  $outpath = substr($outpath,0,-3);
}

open(HANDLE_OUTFILE,">" . $outpath) || die("FATAL ERROR: coluld not create file $outpath");

$currentseq = "";
$currenttag = "";
$currentseqlen = 0;
$dropped = 0; 
$count_seq=0;

while(<HANDLE_INFILE>){

	$line = $_;
	chomp($line);
	
       if(substr($line,0,1) eq ">"){
         $count_seq++;
         if($count_seq>1){
           if($currentseqlen >= $treshold){
             print HANDLE_OUTFILE $currenttag . "\n";
             print HANDLE_OUTFILE $currentseq ;             
           }
           else{
             $dropped++;
           }
         }
         $currenttag = $line;
         $currentseq = "";
         $currentseqlen = 0;
       }
       else{
          $currentseq = $currentseq . $line . "\n" ;
          $currentseqlen = $currentseqlen +length($line);
       }
}

if($currentseqlen >= $treshold){
   print HANDLE_OUTFILE $currenttag . "\n";
   print HANDLE_OUTFILE $currentseq ;             
}
else{
     $dropped++;
}        
           
close(HANDLE_INFILE);
close(HANDLE_OUTFILE);


if( lc(substr($path,-3)) eq ".gz"){

  system("gzip --force $outpath");
  $outpath = $outpath . ".gz";
}

if($overwrite eq "overwrite"){

  system("mv $outpath $path");
}

print "$filecount:processed $path\tcontigs=$count_seq\tdropped=$dropped\n"; 

}


######### SUBS ####################################################################################


sub addtagtopath(){
  my $path = shift(@_);
  my $tag = shift(@_);

  my $newpath=$path . $tag;
  my $extension = "";


  $extension = ".fasta.gz";
  if(substr(lc($path),-1*length($extension)) eq $extension   ){
     $newpath = substr($path,0,-1*length($extension)) . $tag . substr($path,-1*length($extension));
     return $newpath;
  }
  
  $extension = ".fasta";
  if(substr(lc($path),-1*length($extension)) eq $extension   ){
     $newpath = substr($path,0,-1*length($extension)) . $tag . substr($path,-1*length($extension));
    
     return $newpath;
  }
  $extension = ".fa.gz";
  if(substr(lc($path),-1*length($extension)) eq $extension   ){
     $newpath = substr($path,0,-1*length($extension)) . $tag . substr($path,-1*length($extension));
     return $newpath;
  } 
  $extension = ".fa";
  if(substr(lc($path),-1*length($extension)) eq $extension   ){
     $newpath = substr($path,0,-1*length($extension)) . $tag . substr($path,-1*length($extension));
     return $newpath;
  }
  
    $extension = ".fna.gz";
  if(substr(lc($path),-1*length($extension)) eq $extension   ){
     $newpath = substr($path,0,-1*length($extension)) . $tag . substr($path,-1*length($extension));
     return $newpath;
  } 
  $extension = ".fna";
  if(substr(lc($path),-1*length($extension)) eq $extension   ){
     $newpath = substr($path,0,-1*length($extension)) . $tag . substr($path,-1*length($extension));
     return $newpath;
  }
  
  
  return $newpath;
}

