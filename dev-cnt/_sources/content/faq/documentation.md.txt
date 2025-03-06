# Contribute to TSMP documentation.

The quickest and easiest way to contribute to this documentation is to click on 
the `Edit on GitHub` button in the top right corner of each subpage. Clicking 
this button will forward you to the specific document (`.md` or `.rst` file) on 
GitHub from which the subpage is rendered.   
This button gives you two options for editing the documentation.

## Option 1  
You just want to make a small change to this particular subpage. To do this, 
click on the pencil icon (<svg aria-hidden="true" focusable="false" role="img" class="octicon octicon-pencil" viewBox="0 0 16 16" width="16" height="1    6" fill="currentColor" style="display: inline-block; user-select: none; vertical-align: text-bottom; overflow: visible;"><path d="M11.013 1.427a1.75 1.75 0 0 1 2.474 0l1.086 1.086a1.75 1.75 0 0 1 0 2    .474l-8.61 8.61c-.21.21-.47.364-.756.445l-3.251.93a.75.75 0 0 1-.927-.928l.929-3.25c.081-.286.235-.547.445-.758l8.61-8.61Zm.176 4.823L9.75 4.81l-6.286 6.287a.253.253 0 0 0-.064.108l-.558 1.953 1.953-    .558a.253.253 0 0 0 .108-.064Zm1.238-3.763a.25.25 0 0 0-.354 0L10.811 3.75l1.439 1.44 1.263-1.263a.25.25 0 0 0 0-.354Z"></path></svg>) 
in the top right. This will be take you to an editor where you can make your 
change. When you are finished, click on `commit changes`, write a meaningful 
commit message, and commit your change.   
Note, that the `master` branch is protected, and you cannot commit to the 
`master` branch directly. Committing your change will automatically create a 
new branch or fork (depending on your permissions in the repository) and open a 
PR. This is necessary to keep the `master` branch clean, and to ensure that only 
proven content is added to the manual.

## Option 2
You want to edit several parts of the documentation. To do this, click on the 
small down arrow next to the pencil icon (<svg aria-hidden="true" focusable="false" role="img" class="octicon octicon-pencil" viewBox="0 0 16 16" width="16" height="1    6" fill="currentColor" style="display: inline-block; user-select: none; vertical-align: text-bottom; overflow: visible;"><path d="M11.013 1.427a1.75 1.75 0 0 1 2.474 0l1.086 1.086a1.75 1.75 0 0 1 0 2    .474l-8.61 8.61c-.21.21-.47.364-.756.445l-3.251.93a.75.75 0 0 1-.927-.928l.929-3.25c.081-.286.235-.547.445-.758l8.61-8.61Zm.176 4.823L9.75 4.81l-6.286 6.287a.253.253 0 0 0-.064.108l-.558 1.953 1.953-    .558a.253.253 0 0 0 .108-.064Zm1.238-3.763a.25.25 0 0 0-.354 0L10.811 3.75l1.439 1.44 1.263-1.263a.25.25 0 0 0 0-.354Z"></path></svg>). 
This will open a drop down menu, where you can choose between `Edit in place` 
and `github.dev`. `Edit in place` is the default and is used in option 1 above. 
To edit several parts of the documentation, click on `github.dev`. This will 
open the github-dev environment, which is basically a web-based Visual Studio 
Code (VSC) code editor.   

In the left navigation bar you can brows all files of the repository and apply 
your changes. For each file you aply changes at, you will see a increasing 
counter at the left navigation bar with the branch icon. To apply those changes 
click on this branch icon at the left navigation bar. A list of all your changes 
will be shown you, where you can hover each individual change, and choose to 
discard ar stage the change.     
When you are finished, write a meaningful commit message in the text field above 
`commit and push`  and finalize your change by clicking `commit and push`.    
Note that the master branch is protected, and depending on your permissions you 
may not be able to commit to the master branch directly. If this is the case, a 
warning will appear asking you how to proceed. Choose to create a new branch and 
you will be asked to choose a new branch name. Next, Github will ask you to 
switch to the new branch, which you should do if you want to continue editing 
files. When you are done, create a PR for your branch so that your changes can 
be merged into the `master' branch and become part of the documentation. 

## Option 3
As well as web-based options, you can also edit the documentation locally. To do 
this, you need to fork the TSMP repo, create a local clone of the fork, and 
update / change the documentation locally. When you are finished, commit your 
changes with a meaningful commit message, push the new commits to your remote 
fork and create a pull request. 

For your local changes, all you need is a text editor to modify the individual 
`.md` files that make up the documentation. However, to see if your changes 
render as expected, we recommend that you use sphinx to build the entire 
documentation. Do not use a markdown editor for this, as the markdown flavour 
used may differ from that used by sphinx, leading to different results. 

The internal structure of the documentation is quite simple. All files related 
to the documentation are located in `doc/`.   
`conf.py` controls the behaviour of sphinx, `index.rst` is the entry point to 
the documentation, and `content/` contains all the individual `.md`files.  
Building the documentation is done by moving to `doc/` and running the command: 
```
sphinx-build -a . _build
``` 
This will create the documentation in `doc/_build/`. 
Simply browse to this directory and open `index.html` which should show you the 
locally rendered documentation in your default web browser.

## MYST

For referencing we use the extention [myst_parser](https://myst-parser.readthedocs.io/en/latest/).   
To refere to internal pages, you can use the default linking syntax `[label](link)` with relative paths to the internal document. Refering to a specific headline within an internal document you have to use the correct myst ancor. Showing a list of available ancors for a specific document can be dones with:
```
myst-anchors -l auxiliaryscripts.md
>> <h1 id="auxiliary-scripts"></h1>
>> <h2 id="aux_migratefromscratchsh"></h2>
>> <h2 id="aux_untarmanytarssh"></h2>
>> <h2 id="aux_restagetapesh"></h2>
>> <h2 id="aux_gzipsh-and-aux_gunzipsh"></h2>
>> <h2 id="aux_sha512sumsh"></h2>
```
Than use the following link syntax:
```
[aux_UnTarManyTars.sh](REL/PATH/TO/auxiliaryscripts.md#aux_untarmanytarssh)
