# Git instructions overview
You can either clone the repo just to use it (see below in "Clone the repo" section) or pull the master developement branch to push/pull changes (see below in "Initialize Repo" section). If you want to develop your own generation files, in the `generate` sub-directory, feel free to push into master (see below in "Push master" section). If you want to develop a new measurement, molecule, or otherwise modify the files found in this repo, make a new branch (see below in "Branch Repo" section).

## Initialization git:
You need to follow these instructions to do anything with this repo (you can do some stuff on gitlab.com if you have access):
1. Generate an ssh key pair from gitlab.com. see [here](https://docs.gitlab.com/ee/ssh/)
2. Upload the public key to gitlab.com.
3. Add or make an ssh config (`~/.ssh/config`) using the private key (for example `~/.ssh/gitlab_com_rsa`):
```
	#GitLab.com
	Host gitlab.com
		Preferredauthentication publickey
		Identity ~/.ssh/gitlab_com_rsa
```
4. Add the information about git with:
```
	git config --global user.name "Full Name"
	git config --global user.email "user@memphis.edu"
```

## Clone the repo (Using)
If you just want to clone the project to run simulations, follow these instructions:
1. Make sure you initialize git (see "Initialize git" above).
2. Clone the repo:
```
	mkdir mpd-md
	cd mpd-md
	git clone git@gitlab.com:ericspan/mpd-md.git
```

3. Everything, except push and pull usage, is available in this `mpd-md` directory.

## Initialize Repo (Developement)
In order to pull the repo for developement purposes, follow these instructions:
1. Make sure you initialize git (see "Initialize git" above).
2. Pull the repo:
```
	mkdir mpd-md
	cd mpd-md
	git init
	git remote add origin git@gitlab.com:ericspan/mpd-md.git
	git pull origin master
```

3. Everything is in this `mpd-md` directory.

## Push master
If you have a change that you want to push into the main repo, and are absolutely sure you need to, follow these instructions:
1. Make sure you initialize git (see "Initialize git" above).
2. Pull the repo (see "Initialize repo" above).
3. Add your files to wherever they need to be (top, generate, models, etc...).
4. Test your additions, or at least try to compile them.
5. Run the following to add, commit, and push the changes to master (liposome.cpp for example):
```
	git add generate/liposome.cpp
	git commit -m 'Updated liposome.cpp'
	git push origin master
```
6. If you have multiple changes, but are not sure which files you changed, run `git diff` to list them. You can subsequently run everything in step 5 until all the changes are pushed to the master branch.

## Create new branch
To create a new branch, follow these instructions:
To be decided

## Push to a branch
To push into a different branch, follow these instructions:
To be decided, but these are almost identical to "Push master" instructions above.


