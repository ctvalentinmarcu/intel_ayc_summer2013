intel_ayc_summer2013
====================

Intel Accelerate Your Code Contest Summer 2013

Contest page:

http://intel-software-academic-program.com/pages/contests

Problem description:

http://intel-software-academic-program.com/contests/ayc/2013-summer/problem/Intel_AccelerateYourCode_2013-summer_problem.pdf

==================== Intel Summer 2013 AYC Contest =======================

=== Valentin Marcu === Politechnic University of Bucharest === Romania ===


   I define an edge point (of type R, G, or B) as any pixel p which has
at least one neighbour n whose intensity differs from p by more than a
predefined treshold (maxdif), e.g. abs(p.r - n.r) > maxdif. Since
intensities may slightly vary, the stored edge value is approximated to a
close multiple of maxdif (i.e. (int) intensity / maxdif * maxdif).

   For each template, two special pointers (image.rarest and image.biggest)
indicate which group of edge points is the most / least encountered
throughout the image. One of the least encountered edge points (image.rarest[0])
is remembered and called the template's reference.

   The combination of translation, rotation and uniform scaling is implemented
as an affine transform with coefficients a11, a12, a21 and a22.

   An extra ("virtual") translation is performed so that the template's
reference position becomes (0,0). Using this trick, the transform of
(reference.x,reference.y) is (h,w), where (h,w) is the translated candidate 
starting point of the template's transform in the main image. Thus, any points 
that do not match the reference intensities are discarded before computing any
affine coefficients.

  For each (scale, rotation) pair from a set of discrete candidates, three filters
are applied, each with its own predefined maximum of accepted misses: 
 - the comparison between the template.rarest and their transforms;
 - the comparison between the template.biggest and their transforms;
 - the comparison between all the template points and their transforms.

  If a candidate passes all the three filters, the match function returns true, and
h and w are updated in order to reflect the template's (0,0) - transform.
   

