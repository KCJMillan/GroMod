timer
=====

High Precision Timer Class C++/PHP

C++
===
This class calculates the time of execution of a code snippet in seconds, milliseconds or microseconds.                 It is using the new C++11 standard header file chrono (namespace chrono).

Handled exceptions (where the two conflicting calls to timer are done);

Removed cerr as it was also taking up some time.

getElapsedTime() method gives elapsed time at any instance, whether arbitrarily during the process or final.

PHP
===
I have also added a php script for further use.

just use <?php include 'timer.php' ?>
or <?php require 'timer.php' ?>

USAGE:
$time = new Timer; //declare an instance variable
$time->start(); //call the start function
$time->stop(); //stop function returns the time in seconds
