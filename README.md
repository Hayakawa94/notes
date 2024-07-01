How to set up SSH key

$ ssh-keygen -t ed25519 -C "kadt94@gmail.com"
Generating public/private ed25519 key pair.
Enter file in which to save the key (/c/Users/Khoa.Truong/.ssh/id_ed25519):
/c/Users/Khoa.Truong/.ssh/id_ed25519 already exists.
Overwrite (y/n)? y
Enter passphrase (empty for no passphrase):
Enter same passphrase again:
Passphrases do not match.  Try again.
Enter passphrase (empty for no passphrase):
Enter same passphrase again:
Your identification has been saved in /c/Users/Khoa.Truong/.ssh/id_ed25519
Your public key has been saved in /c/Users/Khoa.Truong/.ssh/id_ed25519.pub
The key fingerprint is:
SHA256:VTRLDMwFwQlGrBwTVi06c7sjYyNpjMOnygomLnLHXZw kadt94@gmail.com
The key's randomart image is:
+--[ED25519 256]--+
|      o==*+O*    |
|     .o.o *o.o   |
|     . = .. .    |
|      * ..       |
|       =So       |
|        E        |
|oo + o . .       |
|O = O * o        |
|*=.* o + .       |
+----[SHA256]-----+

locate to c/Users/Khoa.Truong/.ssh/id_ed25519.pub then past it in 
setting >add new ssh keys


link existing project
make sure terminal is in git bash

echo "# RPMtools" >> README.md
git init
git add README.md
git commit -m "first commit"
git branch -M main
git remote add origin git@github.com:Hayakawa94/RPMtools.git
git push -u origin main
