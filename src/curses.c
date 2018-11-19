/*
 *  curses.c
 *  LIME, The versatile 3D line modeling environment
 *
 *  Created by Christian Brinch on 29/10/08.
 *  Copyright 2006-2017, Christian Brinch,
 *  <brinch@nbi.dk>
 *  Niels Bohr institutet
 *  University of Copenhagen
 *	All rights reserved.
 *
 */

#include "lime.h"
#include <curses.h>
#include <time.h>

void
greetings(){
	initscr();
	printw("*** LIME, The versatile line modeling engine, Ver.1.43 (using %d thread(s))\n*** Copyright 2006--2017, Christian Brinch <brinch@nbi.dk>\n", nthreads);
	refresh();
}

void
screenInfo(){
  move(4,4);  printw("Building grid      :");
  move(4,51); printw("|");
  move(5,4);  printw("Smoothing grid     :");
  move(5,51); printw("|");
  move(7,4);  printw("Statistics         :");
  move(9,4);  printw("Iterations         :");
  move(10,4); printw("Propagating photons:");
  move(10,51);printw("|");
  move(13,4); printw("Ray-tracing model  :");
  move(13,51);printw("|");
  move(4,58); printw("|    Molecular data");
  move(5,58); printw("|");
  move(6,58); printw("|");
  move(7,58); printw("|");
  move(8,58); printw("|");
  move(9,58); printw("|");
  move(10,58); printw("|");
  move(11,58); printw("|");
  move(12,58); printw("|");
  move(13,58); printw("|");
  move(14,58); printw("|");
  refresh();
}

void
done(int line){
	move(line,52); printw(" [ok]");
    refresh();
}

void
progressbar(double percent, int line){
  int i;
  for(i=0;i<(int)(percent*25.);i++){
    move(line,25+i);
    printw("#");
  }
  refresh();
}

void
progressbar2(int prog, double percent, double minsnr, double median){
  move(7,38); printw("                    ");
  move(8,38); printw("                    ");
  if(minsnr<1000){
    move(7,25); printw("Min(SNR)    %3.3f", minsnr);
  } else {
    move(7,25); printw("Min(SNR)    %.3e", minsnr);
  }
  if(median<1000){
    move(8,25);	printw("Median(SNR) %3.3f", median);
  } else {
    move(8,25); printw("Median(SNR) %.3e", median);
  }
  move(9,25+prog); printw("#");
  if(percent<100) {
    move(10,25);	 printw("                         ");
  }
  refresh();
}

void
goodnight(int initime, char filename[80]){
	int runtime=time(0)-initime;
	move(14,4); printw("Output written to %s", filename);
	move(22,0); printw("*** Program ended succesfully               ");
	move(22,58); printw("runtime: %3dh %2dm %2ds", runtime/3600, runtime/60%60, runtime%60);
	move(23,0); printw("*** [Press any key to quit]");
    refresh();
	getch();
	endwin();
}

void
quotemass(double mass){
  move(21,6); printw("Total mass contained in model: %3.2e solar masses", mass);
  refresh();
}



void
warning(char message[80]){
	move(22,0); printw("*** %s\n",message);
	refresh();
}

void
bail_out(char message[80]){
	move(22,0); printw("*** %s",message);
	move(23,0); printw("*** [Press any key to quit]");
    refresh();
	getch();
	endwin();
}

void
collpartmesg(char molecule[90], int partners){
  move(6,61); printw("%.25s", molecule);
  move(7,61); printw("%d coll. partner(s):", partners);

  refresh();
}

void
collpartmesg2(char name[10], int partner){
  move(8,61); printw("%s ",name);
  refresh();
}

void
collpartmesg3(int number, int flag){
  move(10,61); printw("Model provides:");
  move(11,61); printw("%d dens. profile(s)", number);
  if(flag==1) {
    move(13,61); printw("*** Warning! ***");
    move(14,61); printw("Too few density profiles");
  }
  refresh();
}
