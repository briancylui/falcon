# Falcon
A very fast hadronization and detector simulator based on ideas pioneered by Bruce Knuteson, in particular, the use of a lookup table to map events at the parton shower level to events at the reconstruction or analysis levels.

---

Text can be **bold**, _italic_, or ~~strikethrough~~.

## [](#header-2)Acknowledgement

The codes for the Falcon package were originally developed by:
* Sergei Gleyzer, University of Florida and CERN;
* Harrison Prosper, Florida State University; and 
* Omar Zapata, UdeA and ITM.

The official Falcon package can be found [here](http://oproject.org/falcon).  The simulator package is mentioned in an academic publication: [Les Houches 2015: Physics at TeV colliders - new physics working group report](http://inspirehep.net/record/1456803#).

## [](#header-2)Background

The current version of Falcon started out as a simulator program that implements a lookup table to record exact matches from parton-level showers to particle-level jets based on geometric proximity, essentially doing what is idiomatically called jet reconstruction.  However, there are two favorable elements missing:
* Compatibility with TMVA
* Learning jet reconstruction rules

## [](#header-2)Project Objectives

In order to improve Falcon on these two areas, this Google Summer of Code 2017 project aims to deploy existing TMVA tools to empower Falcon with the ability to learn the mapping from parton-level showers to particle-level jets.  With this learning capability, TMVA engineers can then scale up the Falcon package to take care of more jet reconstruction tasks with different detector configurations.
 
## [](#header-2)Implementations

Among a basket of plausible methods, implementing a neural network with the MLP (Multi-Layer Perceptron) method of TMVA seems to stand out with its flexibility and fast learning ability.  Users can customize their own settings to train their neural networks for optimal jet reconstruction.  The program exports an Excel (xml) file that stores the learnt weights from the neural network and offers an option for users to launch the TMVA GUI to display regression metrics with a complete set of generated histograms and scatter plots.

## [](#header-2)Future Plans


Falcon

> This is a blockquote following a header.
>
> When something is important enough, you do it even if the odds are not in your favor.

### [](#header-3)Header 3

```js
// Javascript code with syntax highlighting.
var fun = function lang(l) {
  dateformat.i18n = require('./lang/' + l)
  return true;
}
```

```ruby
# Ruby code with syntax highlighting
GitHubPages::Dependencies.gems.each do |gem, version|
  s.add_dependency(gem, "= #{version}")
end
```

#### [](#header-4)Header 4

*   This is an unordered list following a header.
*   This is an unordered list following a header.
*   This is an unordered list following a header.

##### [](#header-5)Header 5

1.  This is an ordered list following a header.
2.  This is an ordered list following a header.
3.  This is an ordered list following a header.

###### [](#header-6)Header 6

| head1        | head two          | three |
|:-------------|:------------------|:------|
| ok           | good swedish fish | nice  |
| out of stock | good and plenty   | nice  |
| ok           | good `oreos`      | hmm   |
| ok           | good `zoute` drop | yumm  |

### There's a horizontal rule below this.

* * *

### Here is an unordered list:

*   Item foo
*   Item bar
*   Item baz
*   Item zip

### And an ordered list:

1.  Item one
1.  Item two
1.  Item three
1.  Item four

### And a nested list:

- level 1 item
  - level 2 item
  - level 2 item
    - level 3 item
    - level 3 item
- level 1 item
  - level 2 item
  - level 2 item
  - level 2 item
- level 1 item
  - level 2 item
  - level 2 item
- level 1 item

### Small image

![](https://assets-cdn.github.com/images/icons/emoji/octocat.png)

### Large image

![](https://guides.github.com/activities/hello-world/branching.png)


### Definition lists can be used with HTML syntax.

<dl>
<dt>Name</dt>
<dd>Godzilla</dd>
<dt>Born</dt>
<dd>1952</dd>
<dt>Birthplace</dt>
<dd>Japan</dd>
<dt>Color</dt>
<dd>Green</dd>
</dl>

```
Long, single-line code blocks should not wrap. They should horizontally scroll if they are too long. This line should be long enough to demonstrate this.
```

```
The final element.
```
