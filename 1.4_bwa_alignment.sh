#!/bin/bash
import subprocess

mlist = [i for i in glob.glob(mdir) if i != mref]
plist = [i for i in glob.glob(pdir) if i != pref]

subprocess.call("bwa mem " + mref + " " + " ".join(mlist) + " > m_alignment.sam")
subprocess.call("bwa mem " + pref + " " + " ".join(plist) + " > p_alignment.sam")
