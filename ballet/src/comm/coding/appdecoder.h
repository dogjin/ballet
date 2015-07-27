/*
 * Copyright (c) 2014 AUTHORS
 *
 * This program is free software; permission to use, copy, modify and
 * distribute this software and its documentation under the terms
 * of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version is hereby granted.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 */

/*!
  \file appdecoder.h
  \brief Definitions of the APPDecoder class
  \author BTS
*/

#ifndef __ballet_comm_appdecoder_h
#define __ballet_comm_appdecoder_h

// ballet includes
#include "ballet/object.h"
#include "convolutionalcoding.h"

// external includes
#include <armadillo>

namespace ballet
{

    class APPDecoderPrivate;

    /*!\addtogroup comm */
    //!@{

    /*!
      \brief A Posteriori Probability Convolutional Decoder

      \par Syntax
      APPDecoder h \n
      APPDecoder h(trellis) \n
      h.decode(Lc,Lu)

      \par Description
      This class implements a posteriori probability (APP) log-map
      decoding of a convolutional code. 

      \par Example
      \code
      //Trellis description
      Trellis t2;
      t2.numInputSymbols = 2;
      t2.numOutputSymbols = 8;
      t2.numStates = 8;
      t2.nextStates = "0 4;4 0;5 1;1 5;2 6;6 2;7 3;3 7";
      t2.codeOutputs = "0 7;0 7;2 5;2 5;3 4;3 4;1 6;1 6";

      fvec Lc = to_fvec(randn(33));

      fvec Lu(ys.length());
      Lu.zeros();

      //Instantiate object
      APPDecoder h;

      //Set object properties
      h.TrellisStructure = t2;
      h.Algorithm = "Max*";

      fvec lud = logmap.decode(Lu,Lc);
      \endcode
    */
    class APPDecoder : public balletObject
    {

        /// @cond INTERNAL
        BALLET_DECLARE_PRIVATE(APPDecoder)
        /// @endcond

    public:

        //! Default constructor
        APPDecoder();

        //! Constructor with trellis
        APPDecoder(const Trellis &trellis);

        /*!
          \brief Decode convolutional code using the a posteriori probability method

          \param[in] Lu         Vector of a priori log-likelihood input bit probabilities
          \param[in] Lc         Vector of a priori log-likelihood coded bit probabilities
          \param[out] lcd       Vector of a posteriori log-likelihood coded bit probabilities

          \par Description
          decode performs APP decoding. The input \b Lu is the sequence of log-likelihoods
          of encoder input data bits. The input \b Lc is the sequence of log-likelihoods
          of encoded bits. Negative values are considered to be zeros and positive values
          are considered to be ones. The method returns \b lud and \b lcd, the updated version 
          of \b Lu, and \b Lc respectively, based on information about the encoder (i.e. trellis,
          termination method, etc.)

        */
        arma::mat decode(const arma::mat &Lu, const arma::mat &Lc, arma::mat *lcd = 0);

    public:

        /*!
          \brief Convolutional code trellis

          The trellis description of the convolutional code.
          The default is the result of poly2trellis(7,[171 133]).
        */
        Trellis TrellisStructure; 

        /*!
          \brief Decoding algorithm (True APP, Max*, Max)

          Specify the decoding algorithm that the object uses
          as one of (True APP, Max*, Max). The default is Max*.
        */
        std::string Algorithm;

        /*!
          \brief Termination method of encoded frame (Truncated, Terminated)

          Specify whether the encoded frame is terminated or truncated.
        */
        std::string TerminationMethod;

        /*!
          \brief Number of scaling bits used to avoid losing precision

          Specify the number of bits the decoder uses to scale the input
          data to avoid losing precision during the computations. The
          default value is 3. The decoder multiplies the input by
          2^NumScalingBits and divides the pre-output by the same factor.
        */
        int NumScalingBits;             ///< Number of bits to scale input to avoid losing precision

    };

    ///@}

};

#endif
