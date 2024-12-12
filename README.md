How to set up git repo
1. Create a new repo in GitHub.com
2. Select SSH 
3. in R create your project with GitHub box tick
4. create a test.R file then do initial comit
5. Copy and paste this in the terminal
		git remote add origin <your SSH>
		git branch -M main
		git push -u origin main

######## SSH KEY is required if this is your first tie set up github
1. paste this in your terminal  ssh-keygen -t ed25519 -C <email address associate with <github account>
2. got to c/Users/Khoa.Truong/.ssh/id_ed25519.pub (this is the default location), but can cd to another location
3. on GitHub.com go to setting>SSH and GPG KEYS then add new SSH Key
4. paste the key  from id_ed25519.pub to the box. which looks like this ssh-ed25519 AAAAC3NzaC1lZDI1NTE5AAAAIL6hZXdfyN2dlIym68LMHJl0o/vfWiQf8ik5lvF1R1Nr <email address associate with github account>
